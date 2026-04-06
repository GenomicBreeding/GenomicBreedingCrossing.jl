"""
    fit_gp_models(genomes::Genomes, phenomes::Phenomes; GB_model::Function = GenomicBreedingModels.bayesa, save::Bool = true, verbose::Bool = true)::Tuple{Dict{String,Fit}, String}

Fit genomic prediction models to phenotypic data for multiple traits.

# Arguments
- `genomes::Genomes`: Genomic data structure containing genotype information.
- `phenomes::Phenomes`: Phenotypic data structure containing trait measurements across individuals.

# Keyword Arguments
- `GB_model::Function`: Genomic prediction model function to fit (default: `GenomicBreedingModels.bayesa`).
- `save::Bool`: Whether to save fitted models to disk (default: `true`).
- `verbose::Bool`: Whether to print progress messages (default: `true`).

# Returns
- `Tuple{Dict{String,Fit}, String}`: A tuple containing:
  - Dictionary mapping trait names to their fitted model objects.
  - Filename of the saved JLD2 file (empty string if `save=false`).

# Details
For each trait in the phenomes, the function:
1. Extracts phenotypic values and removes rows with missing, infinite, or NaN values.
2. Subsets the genomic and phenotypic data to include only valid entries.
3. Fits the specified genomic prediction model.
4. Optionally saves results after each trait is processed.

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> (fits, fname_fits_jld2) = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols, save=false, verbose=false);

julia> (length(fits) == length(phenomes.traits)) && (fname_fits_jld2 == "")
true

julia> length(unique(vcat([v.b_hat_labels for (k, v) in fits]...))) == unique([sum(.!isnan.(v.b_hat)) for (k, v) in fits])[1]
true
```
"""
function fit_gp_models(
    genomes::Genomes,
    phenomes::Phenomes;
    GB_model::Function = GenomicBreedingModels.bayesa,
    save::Bool = true,
    verbose::Bool = true,
)::Tuple{Dict{String,Fit},String}
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false); GB_model::Function = GenomicBreedingModels.ols; save::Bool = true; verbose::Bool = true
    fits = Dict()
    fname_fits_jld2 = save ? "fits-$(string(nameof(GB_model)))-$(Dates.now()).jld2" : ""
    for trait in phenomes.traits
        # trait = phenomes.traits[1]
        if verbose
            println(
                "Fitting \"$(nameof(GB_model))\" model on \"$trait\" and saving as \"$fname_fits_jld2\"",
            )
        end
        G, P = begin
            idx_traits = findall(phenomes.traits .== trait)
            y = phenomes.phenotypes[:, idx_traits][:, 1]
            idx_entries = findall(.!ismissing.(y) .&& .!isinf.(y) .&& .!isnan.(y))
            P = slice(phenomes, idx_entries = idx_entries, idx_traits = idx_traits)
            merge(genomes, P)
        end
        fit = GB_model(genomes = G, phenomes = P, verbose = verbose)
        fits[trait] = fit
        if save
            # Save per trait in case we fail mid-loop
            JLD2.save(fname_fits_jld2, fits)
        end
    end
    (fits, fname_fits_jld2)
end


"""
    predict_gebvs(fit::Fit, genomes::Genomes; idx_entries::Union{Nothing,Vector{Int64}}=nothing)::Phenomes

Predict Genomic Estimated Breeding Values (GEBVs) for a set of genotypes using a fitted model.

# Arguments
- `fit::Fit`: A fitted genomic prediction model.
- `genomes::Genomes`: The genomic data containing genotypes and entry information.
- `idx_entries::Union{Nothing,Vector{Int64}}=nothing`: Optional vector of indices specifying which entries to predict. 
  If `nothing`, predictions are made for all entries in the genomes object.

# Returns
- `Phenomes`: A phenotype object containing the predicted GEBVs for the specified entries, with a single trait.

# Details
The function creates a new `Phenomes` object with predicted values for the selected entries. The predictions 
are computed using the provided fitted model via `GenomicBreedingModels.predict()`. Entry and population 
information are inherited from the input `genomes` object.

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> (fits, fname_fits_jld2) = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols, save=false, verbose=false);

julia> fit = first(sort(fits)).second;

julia> phenomes_predicted = predict_gebvs(fit, genomes);

julia> gebv = phenomes_predicted.phenotypes[:, 1]; gebv = gebv[.!ismissing.(phenomes.phenotypes[:, 1])];

julia> size(phenomes_predicted.phenotypes) == (length(genomes.entries), 1)
true

julia> phenomes_predicted.traits[1] == fit.trait
true

julia> mean((gebv .- fit.y_true).^2) < 0.0001
true
```
"""
function predict_gebvs(
    fit::Fit,
    genomes::Genomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
)::Phenomes
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false); fits, _ = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols); fit = fits["trait_1"]; idx_entries = nothing
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        Int64.(idx_entries)
    end
    phenomes = Phenomes(n = length(idx_entries), t = 1)
    phenomes.traits = [fit.trait]
    phenomes.entries = genomes.entries
    phenomes.populations = genomes.populations
    phenomes.phenotypes[:, 1] = GenomicBreedingModels.predict(
        fit = fit,
        genomes = genomes,
        idx_entries = idx_entries,
    )
    phenomes
end


"""
    panmixis(genomes::Genomes; idx_entries::Union{Nothing,Vector{Int64}} = nothing, 
             population_name::String = "", n::Int64 = 1_000, seed::Int64 = 42, 
             verbose::Bool = false)::Genomes

Simulate random mating (panmixis) of a population to generate progeny genomes.

This function performs random crossing between individuals from a source population,
creating new progeny with allele frequencies that are the average of randomly selected
parent pairs.

# Arguments
- `genomes::Genomes`: Source genomes containing parent individuals and their allele frequencies
- `idx_entries::Union{Nothing,Vector{Int64}}`: Indices of individuals to use as potential parents. 
  If `nothing`, all individuals in `genomes` are used. Default: `nothing`
- `population_name::String`: Name to assign to the progeny population. Default: `""`
- `n::Int64`: Number of progeny to generate. Default: `1_000`
- `seed::Int64`: Random seed for reproducibility. Default: `42`
- `verbose::Bool`: If `true`, display progress bar during simulation. Default: `false`

# Returns
- `Genomes`: A new `Genomes` object containing the simulated progeny with:
  - Allele frequencies computed as the mean of two randomly sampled parents
  - Unique UUIDs for each progeny entry
  - Same loci structure as the source genomes
  - Assigned population name

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> genomes_panmictic = panmixis(genomes, population_name="random_mated");

julia> mean(abs.(mean(genomes_panmictic.allele_frequencies, dims=1) .- mean(genomes.allele_frequencies, dims=1))) < 0.1
true
```
"""
function panmixis(
    genomes::Genomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    population_name::String = "",
    n::Int64 = 1_000,
    seed::Int64 = 42,
    verbose::Bool = false,
)::Genomes
    # genomes=GenomicBreedingCore.simulategenomes(n=119, l=10_000); idx_entries=nothing; population_name::String = ""; n=1_000; seed=42
    Random.seed!(seed)
    rng = Xoshiro(seed)
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        Int64.(idx_entries)
    end
    idx_0 = sample(idx_entries, n, replace = true)
    idx_1 = sample(idx_entries, n, replace = true)
    progenies = Genomes(n = n, p = length(genomes.loci_alleles))
    ids = [uuid4(rng) for i = 1:n]
    progenies.entries = string.("entry-", ids)
    progenies.populations .= population_name
    progenies.loci_alleles = genomes.loci_alleles
    pb = verbose ? ProgressMeter.Progress(n, "Simulating random mating") : nothing
    for (i, (j, k)) in enumerate(zip(idx_0, idx_1))
        # i = 1; j = idx_0[i]; k = idx_1[i]
        progenies.allele_frequencies[i, :] =
            (genomes.allele_frequencies[j, :] .+ genomes.allele_frequencies[k, :]) ./ 2
        verbose ? ProgressMeter.next!(pb) : nothing
    end
    verbose ? ProgressMeter.finish!(pb) : nothing
    @assert checkdims(progenies) # TODO: emit a better more specific error message here
    progenies
end

"""
    simulate_gebv(
        genomes::Genomes,
        phenomes::Phenomes,
        fits::Dict{String, Fit};
        n_cycles::Int64 = 3,
        population_size::Int64 = 1_000,
        selection_intensity::Float64 = 0.1,
        verbose::Bool = true,
        seed::Int64 = 42,
    )::BreedingPopulations

Simulate genomic estimated breeding values (GEBV) through multiple selection cycles using pre-computed genomic prediction model/s.

This function performs iterative selection and breeding cycles for multiple traits. In each cycle,
individuals are selected based on predicted GEBVs, and a new population is generated through panmixis.
Selection intensity determines the proportion of top individuals selected in each cycle.

# Arguments
- `genomes::Genomes`: Initial genome data for the population.
- `phenomes::Phenomes`: Phenotype data containing trait information.
- `fits::Dict{String, Fit}`: Dictionary mapping trait names to fitted genomic prediction models.
- `n_cycles::Int64 = 3`: Number of selection cycles to perform.
- `population_size::Int64 = 1_000`: Size of the population generated after each selection.
- `selection_intensity::Float64 = 0.1`: Proportion of top individuals to select (0.0 to 1.0).
- `verbose::Bool = true`: Whether to display progress meter during simulation.
- `seed::Int64 = 42`: Random seed for reproducibility in panmixis operations.

# Returns
- `BreedingPopulations`: Structure containing selection distributions and phenotypic distributions for all traits and cycles.

# Details
For each trait:
1. Initial GEBV predictions are computed from the input genomes.
2. For each cycle, the top `selection_intensity` proportion of individuals are selected.
3. A new population is generated via panmixis from selected individuals.
4. GEBVs are re-predicted for the new population.
5. Phenotypic distributions are tracked at each stage.

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> (fits, fname_fits_jld2) = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols, save=false, verbose=false);

julia> bp_0 = simulate_gebv(genomes, phenomes, fits, verbose=false);

julia> bp_1 = simulate_gebv(genomes, phenomes, fits, selection_intensity=0.05, verbose=false);

julia> bp_2 = simulate_gebv(genomes, phenomes, fits, selection_intensity=0.95, verbose=false);

julia> mean(bp_1.distributions[end]) >= mean(bp_0.distributions[end]) > mean(bp_2.distributions[end])
true

julia> var(bp_1.distributions[end]) < var(bp_0.distributions[end]) < var(bp_2.distributions[end])
true
```
"""
function simulate_gebv(
    genomes::Genomes,
    phenomes::Phenomes,
    fits::Dict{String,Fit};
    n_cycles::Int64 = 3,
    population_size::Int64 = 1_000,
    selection_intensity::Float64 = 0.1,
    verbose::Bool = true,
    seed::Int64 = 42,
)::BreedingPopulations
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false); fits, _fname_fits_jld2 = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols); n_cycles::Int64 = 3; population_size::Int64 = 1_000; selection_intensity::Float64 = 0.1; save::Bool = true; verbose::Bool = true; seed::Int64 = 42
    bp, _ = BreedingPopulations(n_cycles = n_cycles, n_traits = length(phenomes.traits))
    pb =
        verbose ?
        ProgressMeter.Progress(
            length(phenomes.traits)*n_cycles,
            "Simulating intra-population improvement",
        ) : nothing
    for (j, trait) in enumerate(phenomes.traits)
        # j = 1; trait = phenomes.traits[j]
        fit = fits[trait]
        G = clone(genomes)
        P = predict_gebvs(fit, G)
        bp.traits[j] = trait
        bp.distributions[1, j] = Normal(mean(P.phenotypes[:, 1]), std(P.phenotypes[:, 1])) # cycles x traits
        for t = 1:n_cycles
            # t = 1
            idx_entries = let
                idx = sortperm(P.phenotypes[:, 1], rev = true)
                n = length(idx)
                s = Int64(ceil(selection_intensity * n))
                idx[1:s]
            end
            bp.selections[t, j] = Normal(
                mean(P.phenotypes[idx_entries, 1]),
                std(P.phenotypes[idx_entries, 1]),
            )
            selected_G = slice(G, idx_entries = idx_entries)
            G = panmixis(
                selected_G,
                n = population_size,
                population_name = "cycle_$t",
                seed = seed,
            )
            P = predict_gebvs(fit, G)
            bp.distributions[t+1, j] =
                Normal(mean(P.phenotypes[:, 1]), std(P.phenotypes[:, 1]))
            verbose ? ProgressMeter.next!(pb) : nothing
        end
        # If we were select on the final cycle just for completeness and balancedness of the BreedingPopulations struct
        idx_entries = let
            idx = sortperm(P.phenotypes[:, 1], rev = true)
            n = length(idx)
            s = Int64(ceil(selection_intensity * n))
            idx[1:s]
        end
        bp.selections[n_cycles+1, j] =
            Normal(mean(P.phenotypes[idx_entries, 1]), std(P.phenotypes[idx_entries, 1]))
    end
    if verbose
        ProgressMeter.finish!(pb)
        plot_simple(bp)
    end
    check_populations(bp)
    bp
end

"""
    simulate_gebv(
        genomes::Genomes, 
        phenomes::Phenomes; 
        GB_model::Function = GenomicBreedingModels.bayesa, 
        n_cycles::Int64 = 3, 
        population_size::Int64 = 1_000, 
        selection_intensity::Float64 = 0.1, 
        save::Bool = true, 
        verbose::Bool = true, 
        seed::Int64 = 42
    )::BreedingPopulations

Simulate genomic estimated breeding values (GEBV) through multiple selection cycles.

This function performs iterative selection and breeding cycles for multiple traits. In each cycle,
individuals are selected based on predicted GEBVs, and a new population is generated through panmixis.
Selection intensity determines the proportion of top individuals selected in each cycle.

# Arguments
- `genomes::Genomes`: Genomic data for the population.
- `phenomes::Phenomes`: Phenotypic data for the population.
- `GB_model::Function`: Genomic breeding model to use for predictions (default: `GenomicBreedingModels.bayesa`).
- `n_cycles::Int64`: Number of selection cycles to simulate (default: 3).
- `population_size::Int64`: Size of the population in each cycle (default: 1,000).
- `selection_intensity::Float64`: Intensity of selection applied (default: 0.1).
- `save::Bool`: Whether to save results to disk (default: true).
- `verbose::Bool`: Whether to print progress information (default: true).
- `seed::Int64`: Random seed for reproducibility (default: 42).

# Returns
- `BreedingPopulations`: A breeding population object containing the results of the simulation across all cycles.

# Details
This function first fits genomic prediction models using the specified model function, then simulates the breeding population through multiple selection cycles based on the estimated breeding values.

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> bp_0 = simulate_gebv(genomes, phenomes, verbose=false);

julia> bp_1 = simulate_gebv(genomes, phenomes, selection_intensity=0.05, verbose=false);

julia> bp_2 = simulate_gebv(genomes, phenomes, selection_intensity=0.95, verbose=false);

julia> mean(bp_1.distributions[end]) >= mean(bp_0.distributions[end]) > mean(bp_2.distributions[end])
true

julia> var(bp_1.distributions[end]) < var(bp_0.distributions[end]) < var(bp_2.distributions[end])
true
```
"""
function simulate_gebv(
    genomes::Genomes,
    phenomes::Phenomes;
    GB_model::Function = GenomicBreedingModels.bayesa,
    n_cycles::Int64 = 3,
    population_size::Int64 = 1_000,
    selection_intensity::Float64 = 0.1,
    save::Bool = true,
    verbose::Bool = true,
    seed::Int64 = 42,
)::BreedingPopulations
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false); GB_model::Function = GenomicBreedingModels.ols; n_cycles::Int64 = 3; population_size::Int64 = 1_000; selection_intensity::Float64 = 0.1; save::Bool = true; verbose::Bool = true; seed::Int64 = 42
    fits, fname_fits_jld2 = fit_gp_models(
        genomes,
        phenomes,
        GB_model = GB_model,
        save = save,
        verbose = verbose,
    )
    simulate_gebv(
        genomes,
        phenomes,
        fits;
        n_cycles = n_cycles,
        population_size = population_size,
        selection_intensity = selection_intensity,
        verbose = verbose,
        seed = seed,
    )
end
