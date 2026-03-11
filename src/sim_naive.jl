

function simulate_naive()
    i = 0.1 # i.e. top 10%
    ρ = 0.5 # selection efficiency, i.e. how good are we at phenotyping and/or predicting GEBVs in case of genomic selection
    n_plants_per_cycle = repeat([1_000], 3) # Number of plants per cycle after cycle 0
    h²_realised = 0.5 # h²_realised = 0.5 # h²_realised = R/S = (mean(y_1) - mean(y_0)) / (mean(y_s) - mean(y_0)) -> mean(y_1) = (h²_realised * (mean(y_s) - mean(y_0))) + mean(y_0)
    y::Vector{Union{Missing,Float64}} = rand(1_000)
    selection_direction_to_the_right::Bool = true

    n_cycles = length(n_plants_per_cycle)
    populations = Dict()
    for t in 1:n_cycles
        # t = 1
        n = n_plants_per_cycle[t]
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
        μ = (h²_realised * ρ * (mean(D_s) - mean(D_0))) + mean(D_0) # mean is a function of the heritability, selection efficiency, selections and base population
        # TODO: comment/explain the intuition and that probably something better may exist
        σ = sqrt(var(D_s) + (var(D_0) * h²_realised))
        D_1 = Normal(μ, σ)
        # For now we fix heritability:
        # # TODO: a better way to update h²_realised? --> random walk? Increase? Decrease?
        # r = 0.25*h²_realised * (mean(D_1) - mean(D_0)) / (mean(D_s) - mean(D_0))
        # push!(h²_realised[trait], h²_realised + r*(1-h²_realised))
        # Plot
        y_0 = rand(D_0, n) |> sort
        y_s = rand(D_s, n) |> sort
        y_1 = rand(D_1, n) |> sort
        d_0 = pdf(D_0, y_0)
        d_0 = (d_0 .- minimum(d_0)) ./ (maximum(d_0) .- minimum(d_0))
        d_s = pdf(D_s, y_s)
        d_s = (d_s .- minimum(d_s)) ./ (maximum(d_s) .- minimum(d_s))
        d_s[1] = 0.0
        d_1 = pdf(D_1, y_1)
        d_1 = (d_1 .- minimum(d_1)) ./ (maximum(d_1) .- minimum(d_1))
        p = UnicodePlots.lineplot(y_0, d_0, title="h²=$(round(h²_realised,digits=2))\nR=($(round(mean(D_1),digits=2)) - $(round(mean(D_0),digits=2))) = $(round(mean(D_1)-mean(D_0),digits=2))", xlim=(minimum(vcat(y_0, y_s, y_1)), maximum(vcat(y_0, y_s, y_1))))
        UnicodePlots.lineplot!(p, y_s, d_s)
        UnicodePlots.lineplot!(p, y_1, d_1)
        display(p)
        # Collect information
        populations[t] = Dict(
            "D_0" => [D_0],
            "D_s" => [D_s],
            "D_1" => [D_1],
        )
    end
    populations
end



# # Plot distributions per cycle per trait
# for trait in traits
#     # trait = traits[1]
#     println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
#     println(trait)
#     display(vcat(populations[trait]["D_1"][1], populations[trait]["D_1"]))
#     # For now we fix heritability:
#     # UnicodePlots.barplot(string.(vcat(0, populations[trait]["cycles"])), h²_realised[trait], title="h² ($trait)")
#     # Plot genetic gain over simulated selection cycles
#     ŷ = rand(populations[trait]["D_0"][1], n) |> sort
#     d = pdf(populations[trait]["D_0"][1], ŷ)
#     d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
#     xlim = ((minimum(ŷ)), maximum(vcat(ŷ, rand(populations[trait]["D_1"][end], n))))
#     p = UnicodePlots.lineplot(ŷ, d, title="$trait (h²=$(round(h²_realised,digits=2)))", width=100, xlim=xlim)
#     for t in populations[trait]["cycles"]
#         ŷ = rand(populations[trait]["D_1"][t], n) |> sort
#         d = pdf(populations[trait]["D_1"][t], ŷ)
#         d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
#         UnicodePlots.lineplot!(p, ŷ, d)
#     end
#     display(p)
# end
# # Save naively simulated genetic gains plot
# # Prepare the figure layout
# fig = CairoMakie.Figure(size=(500, 1_500))
# N = length(traits)
# @assert N == 14
# p = 7;
# q = 2;
# axs = Dict()
# for i in 0:(p-1)
#     for j in 0:(q-1)
#         # i = 0; j = 0        
#         idx = (i * q) + j
#         trait = traits[idx+1]
#         ŷ = rand(populations[trait]["D_0"][1], n) |> sort
#         d = pdf(populations[trait]["D_0"][1], ŷ)
#         d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
#         xlim = ((minimum(ŷ)), maximum(vcat(ŷ, rand(populations[trait]["D_1"][end], n))))
#         title = string(split(trait, "|")[2], "\n", replace(replace(replace(split(trait, "|")[3], "NUE" => "highN"), "WUE" => "drought"), "Control" => "lowN"))
#         axs[trait] = CairoMakie.Axis(
#             fig[i, j],
#             title="$title\n(h²=$(round(h²_realised,digits=2)))",
#             limits=(xlim, nothing)
#         )
#         CairoMakie.lines!(axs[trait], ŷ, d)
#     end
# end
# # Plot
# for trait in traits
#     # trait = traits[1]
#     for t in populations[trait]["cycles"]
#         ŷ = rand(populations[trait]["D_1"][t], n) |> sort
#         d = pdf(populations[trait]["D_1"][t], ŷ)
#         d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
#         CairoMakie.lines!(axs[trait], ŷ, d)
#     end
# end
# CairoMakie.save("sample_selection_simulations.svg", fig)

# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
