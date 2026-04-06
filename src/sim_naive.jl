"""
    simulate_naive(phenomes::Phenomes; i::Float64=0.1, ρ::Float64=0.5, 
                   n_plants_per_cycle::Vector{Int64}=repeat([1_000], 5), 
                   h²_realised::Float64=0.5, selection_direction_to_the_right::Bool=true,
                   verbose::Bool=true)::BreedingPopulations

Simulate naive (phenotypic) breeding selection across multiple cycles.

This function models the response to selection in a breeding population using phenotypic selection.
It tracks the distribution of traits across selection cycles, accounting for heritability and 
selection efficiency.

# Arguments
- `phenomes::Phenomes`: The phenotypic data containing trait values for all individuals
- `i::Float64`: Selection intensity (proportion of individuals selected), default 0.1 (top 10%)
- `ρ::Float64`: Selection efficiency (0-1), accounting for accuracy of phenotyping/prediction, default 0.5
- `n_plants_per_cycle::Vector{Int64}`: Number of plants selected per cycle, default 1000 per cycle
- `h²_realised::Float64`: Realized heritability (0-1), default 0.5
- `selection_direction_to_the_right::Bool`: Direction of selection; true for maximization, false for minimization, default true
- `verbose::Bool`: Print progress and plot results, default true

# Returns
- `BreedingPopulations`: A structure containing distributions and selections for each trait across cycles

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames)
julia> _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true, verbose=false);

julia> bp_0 = simulate_naive(phenomes, n_plants_per_cycle=repeat([1_000], 5), verbose=false);

julia> bp_1 = simulate_naive(phenomes, n_plants_per_cycle=repeat([1_000], 5), ρ=1.00, h²_realised=1.00, verbose=false);

julia> bp_2 = simulate_naive(phenomes, n_plants_per_cycle=repeat([1_000], 5), ρ=0.50, h²_realised=0.20, verbose=false);

julia> mean(bp_1.distributions[end]) > mean(bp_0.distributions[end]) > mean(bp_2.distributions[end])
true
```
"""
function simulate_naive(
    phenomes::Phenomes;
    i::Float64 = 0.1, # selection intensity
    ρ::Float64 = 0.5, # selection efficiency, i.e. how good are we at phenotyping and/or predicting GEBVs in case of genomic selection
    n_plants_per_cycle::Vector{Int64} = repeat([1_000], 5), # Number of plants per cycle after cycle 0
    h²_realised::Float64 = 0.5, # h²_realised = 0.5 # h²_realised = R/S = (mean(y_1) - mean(y_0)) / (mean(y_s) - mean(y_0)) -> mean(y_1) = (h²_realised * (mean(y_s) - mean(y_0))) + mean(y_0)
    selection_direction_to_the_right::Bool = true,
    verbose::Bool = true,
)::BreedingPopulations
    # _, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); i = 0.1; ρ = 0.5; n_plants_per_cycle = repeat([1_000], 5); h²_realised = 0.5; selection_direction_to_the_right::Bool = true; verbose::Bool = true
    n_cycles = length(n_plants_per_cycle)
    n_traits = length(phenomes.traits)
    # Initialise with cycle 0 distributions
    bp, _ = BreedingPopulations(n_cycles = n_cycles, n_traits = n_traits)
    for (j, trait) in enumerate(phenomes.traits)
        # j = 1; trait = phenomes.traits[j]
        y = phenomes.phenotypes[:, phenomes.traits .== trait]
        D_0 = let
            ϕ = filter(x -> !ismissing(x) && !isinf(x) && !isnan(x), y) |> x -> Float64.(x)
            sort!(ϕ)
            Normal(mean(ϕ), std(ϕ))
        end
        lower, upper = if selection_direction_to_the_right
            (percentile(D_0, 100 * (1 - i)), Inf)
        else
            (-Inf, percentile(D_0, 100 * i))
        end
        D_s = truncated(D_0, lower, upper)
        bp.traits[j] = trait
        bp.distributions[1, j] = D_0
        bp.selections[1, j] = D_s
    end
    # Across cycles
    if verbose
        pb = ProgressMeter.Progress(n_cycles, "Simulating selection (naive)")
    end
    for t = 1:n_cycles
        # t = 1
        n = n_plants_per_cycle[t]
        bp.cycles[t+1] = t
        for (j, trait) in enumerate(phenomes.traits)
            # j = 1; trait = phenomes.traits[j]
            D_0_prev = bp.distributions[t, j]
            D_s_prev = bp.selections[t, j]
            μ = (h²_realised * ρ * (mean(D_s_prev) - mean(D_0_prev))) + mean(D_0_prev) # mean is a function of the heritability, selection efficiency, selections and base population
            σ = sqrt(var(D_s_prev) + (var(D_0_prev) * h²_realised)) # most probably something better may exist
            D_0 = Normal(μ, σ)
            lower, upper = if selection_direction_to_the_right
                (percentile(D_0, 100 * (1 - i)), Inf)
            else
                (-Inf, percentile(D_0, 100 * i))
            end
            D_s = truncated(D_0, lower, upper)
            bp.distributions[t+1, j] = D_0
            bp.selections[t+1, j] = D_s
            if verbose
                ProgressMeter.next!(pb)
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        plot_simple(bp)
    end
    check_populations(bp)
    bp
end
