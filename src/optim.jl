"""
    sim(;n::Int64=500, l::Int64=10_000)::Tuple{Genomes, Phenomes}

Simulate genomic data and phenotypic observations for a breeding population.

# Arguments
- `n::Int64=500`: Number of individuals in the population.
- `l::Int64=10_000`: Number of genetic loci (markers/SNPs) to simulate.

# Returns
- `Tuple{Genomes, Phenomes}`: A merged dictionary containing both genomic information and phenotypic data for all simulated individuals.

# Description
This function performs an integrated genomic simulation workflow:
1. Generates random genomes for `n` individuals with `l` loci
2. Simulates phenotypic trials with additive, dominance, and epistatic effects
3. Extracts phenotypic measurements from trial results
4. Merges genomic and phenotypic data into a single data structure

# Details
The trial simulation uses a fixed genetic architecture with effects distributed as:
- Additive effects: 50%
- Dominance effects: 40%
- Epistatic effects: 10%

Trials are simulated across a single year, season, harvest, site, and replication.
See GenomicBreedingCore.jl documentation if you wish more fine-tuned simulations.

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes, phenomes = sim(n=10, l=1_000);

julia> size(genomes.allele_frequencies)
(10, 1000)

julia> size(phenomes.phenotypes)
(10, 1)
```

"""
function sim(; n::Int64 = 500, l::Int64 = 10_000)::Tuple{Genomes,Phenomes}
    # n = 500; l = 10_000
    genomes = GenomicBreedingCore.simulategenomes(n = n, l = l, verbose = false)
    trials, simulated_effects = GenomicBreedingCore.simulatetrials(
        genomes = genomes,
        f_add_dom_epi = [0.5 0.4 0.1;],
        n_years = 1,
        n_seasons = 1,
        n_harvests = 1,
        n_sites = 1,
        n_replications = 1,
        verbose = false,
    )
    phenomes = GenomicBreedingCore.extractphenomes(trials)
    merge(genomes, phenomes)
end

"""
    rank_entries(genomes::Genomes, phenomes::Phenomes; idx_trait::Int64=1)::Vector{Int64}

Rank genome entries based on phenotypic values for a specified trait in descending order.

# Arguments
- `genomes::Genomes`: A Genomes struct containing genomic information for entries.
- `phenomes::Phenomes`: A Phenomes struct containing phenotypic data for entries.
- `idx_trait::Int64=1`: The index of the trait to rank entries by (default: 1).

# Returns
- `Vector{Int64}`: A vector of entry indices sorted by phenotypic values in descending order, excluding entries with missing phenotypic data.

# Throws
- `String`: An error message if the number of entries in `genomes` and `phenomes` do not match.

# Notes
- Only entries with non-missing phenotypic values for the specified trait are included in the ranking.
- Missing values in the phenotype data are automatically filtered out before ranking.

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes, phenomes = sim();

julia> idx = rank_entries(genomes, phenomes, idx_trait=1);

julia> idx == sortperm(phenomes.phenotypes[:, 1], rev=true)
true
```
"""
function rank_entries(
    genomes::Genomes,
    phenomes::Phenomes;
    idx_trait::Int64 = 1,
)::Vector{Int64}
    # genomes, phenomes = sim(); idx_trait = 1; phenomes.phenotypes[1, idx_trait] = missing
    if genomes.entries != phenomes.entries
        throw(
            "Please make sure that the Genomes and Phenomes struct have been merged and so correspond to each other's entries.",
        )
    end
    idx_no_missing = findall(.!ismissing.(phenomes.phenotypes[:, idx_trait]))
    y::Vector{Float64} = phenomes.phenotypes[idx_no_missing, idx_trait]
    idx_sorted = sortperm(y, rev = true)
    idx_no_missing[idx_sorted]
end

"""
    two_way_cross_predictions(
        genomes::Genomes,
        phenomes::Phenomes;
        idx_trait::Int64=1,
        n_parents::Int64=5,
        n::Int64 = 1_000,
        seed::Int64 = 42,
        ϵ::Float64 = 0.001,
        GB_model::Function = GenomicBreedingModels.bayesa,
        verbose::Bool = false,
    )::DataFrame

Evaluate all pairwise crosses between top-ranked parents and predict offspring phenotypes.

This function identifies the top `n_parents` individuals based on their trait values, generates all possible 
pairwise crosses between these parents, and predicts genomic estimated breeding values (GEBVs) for the 
resulting offspring populations.

# Arguments
- `genomes::Genomes`: Genotype data for all individuals
- `phenomes::Phenomes`: Phenotype data for all individuals
- `idx_trait::Int64=1`: Index of the trait to optimize (default: 1)
- `n_parents::Int64=5`: Number of top-ranked parents to select for crossing (default: 5)
- `n::Int64=1_000`: Number of offspring to simulate per cross (default: 1,000)
- `seed::Int64=42`: Random seed for reproducibility (default: 42)
- `ϵ::Float64=0.001`: Recombination rate parameter for paired crosses (default: 0.001)
- `GB_model::Function=GenomicBreedingModels.bayesa`: Genomic prediction model to use (default: BayesA)
- `verbose::Bool=false`: Print progress updates (default: false)

# Returns
- `DataFrame`: Results containing columns:
  - `parent_1::Vector{String}`: Name/ID of first parent
  - `parent_2::Vector{String}`: Name/ID of second parent
  - `ϕ_min::Vector{Float64}`: Minimum GEBV in offspring
  - `ϕ_max::Vector{Float64}`: Maximum GEBV in offspring
  - `ϕ_mean::Vector{Float64}`: Mean GEBV in offspring
  - `ϕ_var::Vector{Float64}`: Variance of GEBVs in offspring

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes, phenomes = sim();

julia> df = two_way_cross_predictions(genomes, phenomes; idx_trait=1, n_parents=5, GB_model = GenomicBreedingModels.ols);

julia> size(df)
(10, 6)
```
"""
function two_way_cross_predictions(
    genomes::Genomes,
    phenomes::Phenomes;
    idx_trait::Int64 = 1,
    n_parents::Int64 = 5,
    n::Int64 = 1_000,
    seed::Int64 = 42,
    ϵ::Float64 = 0.001,
    GB_model::Function = GenomicBreedingModels.bayesa,
    verbose::Bool = false,
)::DataFrame
    # genomes, phenomes = sim(); idx_trait = 1; n_parents::Int64=5; n::Int64 = 1_000; seed::Int64 = 42; ϵ::Float64 = 0.001; GB_model::Function = GenomicBreedingModels.ols; verbose::Bool = true
    Σ = cov(genomes.allele_frequencies)
    idx = rank_entries(genomes, phenomes, idx_trait = idx_trait)[1:n_parents]
    fit = let
        fits, _ = fit_gp_models(
            genomes,
            slice(phenomes, idx_entries = idx, idx_traits = [idx_trait]);
            GB_model = GB_model,
            save = false,
            verbose = verbose,
        )
        fits[string.(keys(fits))[1]]
    end
    m = Int64(round((n_parents^2 - n_parents) / 2))
    parent_1::Vector{String} = repeat([""], m)
    parent_2::Vector{String} = repeat([""], m)
    ϕ_min::Vector{Float64} = repeat([NaN], m)
    ϕ_max::Vector{Float64} = repeat([NaN], m)
    ϕ_mean::Vector{Float64} = repeat([NaN], m)
    ϕ_var::Vector{Float64} = repeat([NaN], m)
    pb =
        verbose ?
        ProgressMeter.Progress(
            m,
            desc = "Two-way cross predictions for $m pairs of the top-performing genotypes",
        ) : nothing
    Threads.@threads for i = 1:(length(idx)-1)
        for j = (i+1):length(idx)
            # i = 1; j = 2
            genomes_1 = slice(genomes, idx_entries = [idx[i]])
            genomes_2 = slice(genomes, idx_entries = [idx[j]])
            Γ::Genomes = paired_cross(
                genomes_1,
                genomes_2;
                Σ = Σ,
                population_name = string(genomes_1.entries[1], "▓", genomes_2.entries[1]),
                n = n,
                seed = seed,
                ϵ = ϵ,
                verbose = false,
            )
            Φ::Phenomes = predict_gebvs(fit, Γ)
            parent_1[i] = genomes_1.entries[1]
            parent_2[i] = genomes_2.entries[1]
            ϕ_min[i] = minimum(Φ.phenotypes[:, 1])
            ϕ_max[i] = maximum(Φ.phenotypes[:, 1])
            ϕ_mean[i] = mean(Φ.phenotypes[:, 1])
            ϕ_var[i] = var(Φ.phenotypes[:, 1])
            verbose ? ProgressMeter.next!(pb) : nothing
        end
    end
    verbose ? ProgressMeter.finish!(pb) : nothing
    DataFrame(
        parent_1 = parent_1,
        parent_2 = parent_2,
        ϕ_min = ϕ_min,
        ϕ_max = ϕ_max,
        ϕ_mean = ϕ_mean,
        ϕ_var = ϕ_var,
    )
end

"""
    synthetic_predictions(genomes::Genomes, phenomes::Phenomes; 
                         idx_trait::Int64=1, n_total_candidates::Int64=10, 
                         n_parents_per_synthetic::Int64=5, n::Int64=1_000, 
                         seed::Int64=42, ϵ::Float64=0.001, 
                         GB_model::Function=GenomicBreedingModels.bayesa, 
                         verbose::Bool=false)::DataFrame

Generate predictions for synthetic crosses by evaluating all combinations of top-performing parental genotypes.

# Arguments
- `genomes::Genomes`: The genomic data containing allele frequencies and genotype information.
- `phenomes::Phenomes`: The phenotypic data corresponding to the genomes.
- `idx_trait::Int64=1`: Index of the trait to evaluate for parent selection.
- `n_total_candidates::Int64=10`: Number of top-performing candidates to consider as potential parents.
- `n_parents_per_synthetic::Int64=5`: Number of parents to combine in each synthetic cross.
- `n::Int64=1_000`: Number of individuals to generate for each synthetic population.
- `seed::Int64=42`: Random seed for reproducibility of synthetic population generation.
- `ϵ::Float64=0.001`: Small tolerance parameter (currently unused).
- `GB_model::Function=GenomicBreedingModels.bayesa`: Genomic breeding model function for GEBV prediction.
- `verbose::Bool=false`: If true, display progress bar and additional information during execution.

# Returns
- `DataFrame`: A table with columns:
  - `parents::Vector{String}`: Concatenated IDs of parent genotypes for each cross.
  - `ϕ_min::Vector{Float64}`: Minimum predicted phenotypic value in the synthetic population.
  - `ϕ_max::Vector{Float64}`: Maximum predicted phenotypic value in the synthetic population.
  - `ϕ_mean::Vector{Float64}`: Mean predicted phenotypic value in the synthetic population.
  - `ϕ_var::Vector{Float64}`: Variance of predicted phenotypic values in the synthetic population.

# Description
This function ranks genotypes based on phenotypic performance, selects the top candidates, 
fits a predictive model, and evaluates all possible combinations of parent crosses using panmixis. 
Genomic estimated breeding values (GEBVs) are predicted for each synthetic population.

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes, phenomes = sim();

julia> df = synthetic_predictions(genomes, phenomes; idx_trait=1, n_total_candidates=3, n_parents_per_synthetic=2, GB_model = GenomicBreedingModels.ols);

julia> size(df)
(3, 5)
```
"""
function synthetic_predictions(
    genomes::Genomes,
    phenomes::Phenomes;
    idx_trait::Int64 = 1,
    n_total_candidates::Int64 = 10,
    n_parents_per_synthetic::Int64 = 5,
    n::Int64 = 1_000,
    seed::Int64 = 42,
    ϵ::Float64 = 0.001,
    GB_model::Function = GenomicBreedingModels.bayesa,
    verbose::Bool = false,
)::DataFrame
    # genomes, phenomes = sim(); idx_trait = 1; n_total_candidates::Int64 = 10; n_parents_per_synthetic::Int64=5; n::Int64 = 1_000; seed::Int64 = 42; ϵ::Float64 = 0.001; GB_model::Function = GenomicBreedingModels.ols; verbose::Bool = true
    idx = rank_entries(genomes, phenomes, idx_trait = idx_trait)[1:n_total_candidates]
    fit = let
        fits, _ = fit_gp_models(
            genomes,
            slice(phenomes, idx_traits = [idx_trait]);
            GB_model = GB_model,
            save = false,
            verbose = verbose,
        )
        fits[string.(keys(fits))[1]]
    end
    idx_parents_combinations = collect(combinations(idx, n_parents_per_synthetic))
    m = length(idx_parents_combinations)
    parents::Vector{String} = repeat([""], m)
    ϕ_min::Vector{Float64} = repeat([NaN], m)
    ϕ_max::Vector{Float64} = repeat([NaN], m)
    ϕ_mean::Vector{Float64} = repeat([NaN], m)
    ϕ_var::Vector{Float64} = repeat([NaN], m)
    pb =
        verbose ?
        ProgressMeter.Progress(
            m,
            desc = "Two-way cross predictions for $m pairs of the top-performing genotypes",
        ) : nothing
    Threads.@threads for i = 1:length(idx_parents_combinations)
        # i = 1; 
        idx_parents = idx_parents_combinations[i]
        id = join(genomes.entries[idx_parents], "▓")
        Γ::Genomes = panmixis(
            genomes,
            idx_entries = idx_parents,
            population_name = id,
            n = n,
            simple_mating_model = false,
            ϵ = ϵ,
            seed = seed,
            verbose = false,
        )
        Φ::Phenomes = predict_gebvs(fit, Γ)
        parents[i] = id
        ϕ_min[i] = minimum(Φ.phenotypes[:, 1])
        ϕ_max[i] = maximum(Φ.phenotypes[:, 1])
        ϕ_mean[i] = mean(Φ.phenotypes[:, 1])
        ϕ_var[i] = var(Φ.phenotypes[:, 1])
        verbose ? ProgressMeter.next!(pb) : nothing
    end
    verbose ? ProgressMeter.finish!(pb) : nothing
    DataFrame(
        parents = parents,
        ϕ_min = ϕ_min,
        ϕ_max = ϕ_max,
        ϕ_mean = ϕ_mean,
        ϕ_var = ϕ_var,
    )
end


"""
    test_cross(genomes::Genomes, phenomes::Phenomes; 
               idx_trait::Int64=1, idx_tester::Int64=1, n::Int64=1_000, 
               seed::Int64=42, ϵ::Float64=0.001, 
               GB_model::Function=GenomicBreedingModels.bayesa, 
               verbose::Bool=false)::DataFrame

Efficiently test crosses between a tester genotype and all other genotypes in a population.

This function performs genomic prediction by crossing a specified tester genotype against each 
other genotype in the population, computing genomic estimated breeding values (GEBVs) for the 
resulting synthetic offspring.

# Arguments
- `genomes::Genomes`: Population of genotypes to evaluate
- `phenomes::Phenomes`: Corresponding phenotypic data for model training
- `idx_trait::Int64=1`: Index of the trait to evaluate
- `idx_tester::Int64=1`: Index of the tester genotype in the population
- `n::Int64=1_000`: Number of synthetic offspring to simulate per cross
- `seed::Int64=42`: Random seed for reproducibility of cross simulations
- `ϵ::Float64=0.001`: Tolerance parameter for cross generation
- `GB_model::Function=GenomicBreedingModels.bayesa`: Genomic prediction model function
- `verbose::Bool=false`: Enable progress reporting

# Returns
- `DataFrame`: Results table with columns:
  - `tester_parents::Vector{String}`: Tester genotype identifier
  - `other_parents::Vector{String}`: Other parent genotype identifier
  - `ϕ_min::Vector{Float64}`: Minimum GEBV of offspring
  - `ϕ_max::Vector{Float64}`: Maximum GEBV of offspring
  - `ϕ_mean::Vector{Float64}`: Mean GEBV of offspring
  - `ϕ_var::Vector{Float64}`: Variance of GEBVs in offspring

# Notes
Parallel processing is used via `Threads.@threads` for computational efficiency.

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes, phenomes = sim(n=20);

julia> df = test_cross(genomes, phenomes, GB_model = GenomicBreedingModels.ols);

julia> size(df)
(19, 6)
```

"""
function test_cross(
    genomes::Genomes,
    phenomes::Phenomes;
    idx_trait::Int64 = 1,
    idx_tester::Int64 = 1,
    n::Int64 = 1_000,
    seed::Int64 = 42,
    ϵ::Float64 = 0.001,
    GB_model::Function = GenomicBreedingModels.bayesa,
    verbose::Bool = false,
)::DataFrame
    # genomes, phenomes = sim(n=20); idx_trait = 1; idx_tester::Int64 = 1; n::Int64 = 1_000; seed::Int64 = 42; ϵ::Float64 = 0.001; GB_model::Function = GenomicBreedingModels.ols; verbose::Bool = true
    Σ = cov(genomes.allele_frequencies)
    fit = let
        fits, _ = fit_gp_models(
            genomes,
            slice(phenomes, idx_traits = [idx_trait]);
            GB_model = GB_model,
            save = false,
            verbose = verbose,
        )
        fits[string.(keys(fits))[1]]
    end
    idx_other_parents = setdiff(collect(1:length(genomes.entries)), [idx_tester])
    m = length(idx_other_parents)
    tester_parents::Vector{String} = repeat([""], m)
    other_parents::Vector{String} = repeat([""], m)
    ϕ_min::Vector{Float64} = repeat([NaN], m)
    ϕ_max::Vector{Float64} = repeat([NaN], m)
    ϕ_mean::Vector{Float64} = repeat([NaN], m)
    ϕ_var::Vector{Float64} = repeat([NaN], m)
    genomes_1 = slice(genomes, idx_entries = [idx_tester])
    pb =
        verbose ?
        ProgressMeter.Progress(
            m,
            desc = "Test-crossing $m genotypes against $(genomes.entries[idx_tester]) tester",
        ) : nothing
    Threads.@threads for i = 1:length(idx_other_parents)
        # i = 1
        idx = idx_other_parents[i]
        Γ::Genomes = paired_cross(
            genomes_1,
            slice(genomes, idx_entries = [idx]);
            Σ = Σ,
            population_name = genomes.entries[idx],
            n = n,
            seed = seed,
            ϵ = ϵ,
            verbose = false,
        )
        Φ::Phenomes = predict_gebvs(fit, Γ)
        tester_parents[i] = genomes.entries[idx_tester]
        other_parents[i] = genomes.entries[idx]
        ϕ_min[i] = minimum(Φ.phenotypes[:, 1])
        ϕ_max[i] = maximum(Φ.phenotypes[:, 1])
        ϕ_mean[i] = mean(Φ.phenotypes[:, 1])
        ϕ_var[i] = var(Φ.phenotypes[:, 1])
        verbose ? ProgressMeter.next!(pb) : nothing
    end
    verbose ? ProgressMeter.finish!(pb) : nothing
    DataFrame(
        tester_parents = tester_parents,
        other_parents = other_parents,
        ϕ_min = ϕ_min,
        ϕ_max = ϕ_max,
        ϕ_mean = ϕ_mean,
        ϕ_var = ϕ_var,
    )
end
