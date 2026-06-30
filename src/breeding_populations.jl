mutable struct BreedingPopulations
    cycles::Vector{Int64}
    traits::Vector{String}
    distributions::Matrix{Distribution} # cycles (including cycle 0) x traits
    selections::Matrix{Distribution} # cycles (including cycle 0) x traits
    function BreedingPopulations(;
        n_cycles::Int64 = 5,
        n_traits::Int64 = 3,
        simulate_genomes_phenomes::Bool = false,
        verbose::Bool = true,
    )::Tuple{BreedingPopulations,Union{Nothing,Tuple{Genomes,Phenomes}}}
        # n_cycles=5; n_traits=3; simulate_genomes_phenomes=false; verbose=false
        cycles::Vector{Int64} = collect(0:n_cycles)
        traits::Vector{String} = ["trait_$i" for i = 1:n_traits]
        distributions::Matrix{Distribution} =
            reshape(repeat([Normal()], (n_cycles + 1) * n_traits), (n_cycles + 1), n_traits)
        selections::Matrix{Distribution} =
            reshape(repeat([Normal()], (n_cycles + 1) * n_traits), (n_cycles + 1), n_traits)
        if !simulate_genomes_phenomes
            return (new(cycles, traits, distributions, selections), nothing)
        else
            # Random unique contents
            genomes = simulategenomes(verbose = verbose)
            trials, _ = simulatetrials(
                genomes = genomes,
                n_years = 1,
                n_seasons = 1,
                n_measurements = 1,
                n_sites = 1,
                n_replications = 1,
                sparsity = 0.10,
                verbose = verbose,
            )
            phenomes = extractphenomes(trials)
            for i = 1:(n_cycles+1)
                for j = 1:n_traits
                    # i = 1; j = 1
                    Random.seed!((i * j) + j)
                    μ = try
                        u = mean(phenomes.phenotypes[:, j]) * i
                        s = std(phenomes.phenotypes[:, j]) * (1/i)
                        rand(Normal(u, s))
                    catch
                        rand(Normal(i))
                    end
                    σ = abs(rand(Normal()))
                    distributions[i, j] = Normal(μ, σ)
                    selections[i, j] = truncated(Normal(μ, σ), 1.75 * μ, Inf)
                end
            end
            return (new(cycles, traits, distributions, selections), (genomes, phenomes))
        end
    end
end

function check_populations(bp::BreedingPopulations)::Nothing
    # bp, _ = BreedingPopulations()
    n_cycles = length(bp.cycles)
    n_traits = length(bp.traits)
    if length(n_cycles) == 0
        throw(ErrorException("No selection cycles found"))
    end
    if length(n_traits) == 0
        throw(ErrorException("No traits found"))
    end
    if length(unique(bp.traits)) < n_traits
        dups = filter(x -> x[2] > 1, countmap(bp.traits))
        throw(ErrorException("Duplicate trait names: $dups"))
    end
    if size(bp.distributions) != (n_cycles, n_traits)
        throw(
            ErrorException(
                "The distributions matrix is incompatible with the number of cycles and/or traits (n_cycles=$n_cycles, n_traits=$n_traits)",
            ),
        )
    end
    if size(bp.selections) != (n_cycles, n_traits)
        throw(
            ErrorException(
                "The selections matrix is incompatible with the number of cycles and/or traits (n_cycles=$n_cycles, n_traits=$n_traits)",
            ),
        )
    end
    nothing
end

function plot_simple(bp::BreedingPopulations; n::Int64 = 10_000)::Nothing
    # bp, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); n = 10_000
    check_populations(bp)
    n_cycles = length(bp.cycles) - 1
    for (j, trait) in enumerate(bp.traits)
        # j = 1; trait = bp.traits[j]
        y_min, y_max = (Inf, -Inf)
        for i = 1:(n_cycles+1)
            ŷ = sort(rand(bp.distributions[i, j], n))
            y_min = minimum(ŷ) < y_min ? minimum(ŷ) : y_min
            y_max = maximum(ŷ) > y_max ? maximum(ŷ) : y_max
        end
        ŷ = sort(rand(bp.distributions[1, j], n))
        f = pdf(bp.distributions[1, j], ŷ)
        f = (f .- minimum(f)) ./ (maximum(f) - minimum(f))
        p = UnicodePlots.scatterplot(
            ŷ,
            f,
            xlim = (y_min, y_max),
            title = trait,
            name = "Cycle 0",
        )
        for i = 2:(n_cycles+1)
            ŷ = sort(rand(bp.distributions[i, j], n))
            f = pdf(bp.distributions[i, j], ŷ)
            f = (f .- minimum(f)) ./ (maximum(f) - minimum(f))
            UnicodePlots.scatterplot!(p, ŷ, f, name = "Cycle $(i-1)")
        end
        display(p)
    end
    nothing
end

function plot_and_save(
    bp::BreedingPopulations;
    n::Int64 = 10_000,
    figure_size::Tuple{Int64,Int64} = (900, 500),
    map_pdf_0_to_1::Bool = true,
    color_scheme::ColorScheme = ColorSchemes.Tam,
    fname_svg::String = "",
)::String
    # bp, (genomes, phenomes) = BreedingPopulations(simulate_genomes_phenomes=true); n::Int64=10_000; figure_size::Tuple{Int64,Int64}=(900, 500); map_pdf_0_to_1::Bool=true; color_scheme::ColorScheme = ColorSchemes.Tam; fname_svg::String = ""
    check_populations(bp)
    n_cycles = length(bp.cycles) - 1
    n_traits = length(bp.traits)
    fname_svg = if fname_svg == ""
        "breeding_populations-$(Dates.now())-$(Int64(ceil(100_000*rand()))).svg"
    else
        fname_svg
    end
    colours = resample(color_scheme, n_cycles + 1, (alpha) -> 0.8)
    fig = CairoMakie.Figure(size = figure_size)
    layout = fig[1, 1] = GridLayout()
    axes = Dict()
    n_rows = Int64(ceil(sqrt(n_traits)))
    n_cols = Int64(ceil(n_traits / n_rows))
    i = 1
    j = 0
    Ys = Dict()
    for (idx_trait, trait) in enumerate(bp.traits)
        # idx_trait = 1; trait = bp.traits[idx_trait]
        i = (j == n_cols) && (i < n_rows) ? i + 1 : i
        j = j < n_cols ? j + 1 : 1
        y_min, y_max = (Inf, -Inf)
        Ys[trait] = Dict()
        for idx_cycle = 1:(n_cycles+1)
            # idx_cycle = 1
            ŷ = sort(rand(bp.distributions[idx_cycle, idx_trait], n))
            y_min = minimum(ŷ) < y_min ? minimum(ŷ) : y_min
            y_max = maximum(ŷ) > y_max ? maximum(ŷ) : y_max
            f = pdf(bp.distributions[idx_cycle, idx_trait], ŷ)
            f = map_pdf_0_to_1 ? (f .- minimum(f)) ./ (maximum(f) - minimum(f)) : f
            Ys[trait]["Cycle $(idx_cycle-1)"] = Dict("ŷ" => ŷ, "f" => f)
        end
        axes[trait] =
            CairoMakie.Axis(layout[i, j], title = trait, limits = ((y_min, y_max), nothing))
    end

    for trait in bp.traits
        # trait = bp.traits[1]
        for idx_cycle = 1:(n_cycles+1)
            x = Ys[trait]["Cycle $(idx_cycle-1)"]["ŷ"]
            y = Ys[trait]["Cycle $(idx_cycle-1)"]["f"]
            CairoMakie.lines!(axes[trait], x, y, color = colours[idx_cycle])
        end
    end
    elems = [
        [MarkerElement(color = x, marker = :circle, markersize = 15, strokecolor = :black)] for x in colours
    ]
    CairoMakie.Legend(fig[1, 2], elems, ["Cycle $t" for t = 0:n_cycles])
    CairoMakie.save(fname_svg, fig)
    fname_svg
end
