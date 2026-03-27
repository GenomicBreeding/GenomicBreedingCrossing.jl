function fit_gp_models(
    genomes::Genomes,
    phenomes::Phenomes;
    GB_model::Function = GenomicBreedingModels.bayesa,
    save::Bool = true,
    verbose::Bool = true,
)::Tuple{Dict{String,Fit}, String}
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); GB_model::Function = GenomicBreedingModels.ols; save::Bool = true; verbose::Bool = true
    fits = Dict()
    fname_fits_jld2 = save ? "fits-$(string(nameof(GB_model)))-$(Dates.now()).jld2" : ""
    for trait in phenomes.traits
        # trait = phenomes.traits[1]
        if verbose
            println("Fitting \"$(nameof(GB_model))\" model on \"$trait\" and saving as \"$fname_fits_jld2\"")
        end
        G, P = begin
            idx_traits = findall(phenomes.traits .== trait)
            y = phenomes.phenotypes[:, idx_traits][:, 1]
            idx_entries = findall(.!ismissing.(y) .&& .!isinf.(y) .&& .!isnan.(y))
            P = slice(phenomes, idx_entries = idx_entries, idx_traits = idx_traits)
            merge(genomes, P)
        end
        fit = GB_model(genomes = G, phenomes = P, verbose=verbose)
        fits[trait] = fit
        if save
            # Save per trait in case we fail mid-loop
            JLD2.save(fname_fits_jld2, fits)
        end
    end
    (fits, fname_fits_jld2)
end


# GEBVs
function predict_gebvs(
    fit::Fit,
    genomes::Genomes;
    idx_entries::Union{Nothing,Vector{Int64}} = nothing,
)::Phenomes
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); fits, _ = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols); fit = fits["trait_1"]; idx_entries = nothing
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        Int64.(idx_entries)
    end
    phenomes = Phenomes(n = length(idx_entries), t = 1)
    phenomes.entries = genomes.entries
    phenomes.populations = genomes.populations
    phenomes.phenotypes[:, 1] = GenomicBreedingModels.predict(
        fit = fit,
        genomes = genomes,
        idx_entries = idx_entries,
    )
    phenomes
end


# Simulate mating
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
    idx_0 = sample(idx_entries, n, replace=true)
    idx_1 = sample(idx_entries, n, replace=true)
    progenies = Genomes(n = n, p = length(genomes.loci_alleles))
    ids = [uuid4(rng) for i in 1:n]
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


# _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true)
# fits, _ = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols)
# bp = simulate_gebv(genomes, phenomes, fits)
# plot_simple(bp)
function simulate_gebv(
    genomes::Genomes,
    phenomes::Phenomes,
    fits::Dict{String, Fit};
    n_cycles::Int64 = 3,
    population_size::Int64 = 1_000,
    selection_intensity::Float64 = 0.1,
    verbose::Bool = true,
    seed::Int64 = 42,
)::BreedingPopulations
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); fits, _fname_fits_jld2 = fit_gp_models(genomes, phenomes, GB_model=GenomicBreedingModels.ols); n_cycles::Int64 = 3; population_size::Int64 = 1_000; selection_intensity::Float64 = 0.1; save::Bool = true; verbose::Bool = true; seed::Int64 = 42
    bp, _ = BreedingPopulations(n_cycles=n_cycles, n_traits=length(phenomes.traits))
    pb = verbose ? ProgressMeter.Progress(length(phenomes.traits)*n_cycles, "Simulating intra-population improvement") : nothing
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
            bp.selections[t, j] = Normal(mean(P.phenotypes[idx_entries, 1]), std(P.phenotypes[idx_entries, 1]))
            selected_G = slice(G, idx_entries = idx_entries)
            G = panmixis(
                selected_G,
                n = population_size,
                population_name = "cycle_$t",
                seed = seed,
            )
            P = predict_gebvs(fit, G)
            bp.distributions[t+1, j] = Normal(mean(P.phenotypes[:, 1]), std(P.phenotypes[:, 1]))
            verbose ? ProgressMeter.next!(pb) : nothing
        end
        # If we were select on the final cycle just for completeness and balancedness of the BreedingPopulations struct
        idx_entries = let
            idx = sortperm(P.phenotypes[:, 1], rev = true)
            n = length(idx)
            s = Int64(ceil(selection_intensity * n))
            idx[1:s]
        end
        bp.selections[n_cycles+1, j] = Normal(mean(P.phenotypes[idx_entries, 1]), std(P.phenotypes[idx_entries, 1]))
    end
    verbose ? ProgressMeter.finish!(pb) : nothing
    check_populations(bp)
    bp
end
# _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true)
# bp = simulate_gebv(genomes, phenomes, GB_model=GenomicBreedingModels.ols)
# plot_simple(bp)
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
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); GB_model::Function = GenomicBreedingModels.ols; n_cycles::Int64 = 3; population_size::Int64 = 1_000; selection_intensity::Float64 = 0.1; save::Bool = true; verbose::Bool = true; seed::Int64 = 42
    fits, fname_fits_jld2 = fit_gp_models(
        genomes,
        phenomes,
        GB_model=GB_model,
        save=save,
        verbose=verbose,
    )
    simulate_gebv(
        genomes,
        phenomes,
        fits;
        n_cycles=n_cycles,
        population_size=population_size,
        selection_intensity=selection_intensity,
        verbose=verbose,
        seed=seed,
    )
end
