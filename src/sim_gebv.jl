
#############################################
### Trying to include genomes simulations ###
#############################################
traits = names(df)[2:end]
G_0 = readjld2(Genomes, fname="/group/pasture/forages/Ryegrass/STR/NUE_WUE_2022-23/NUE_WUE_halfsibs/genotype/genomes_v2.jld2")
# G_0 = slice(G_0, idx_loci_alleles=collect(1:13:length(G_0.loci_alleles)))
G_0.populations .= "cycle_0"
dimensions(G_0)
# fits = try
#     JLD2.load("fits.jld2")
# catch
#     fits = Dict()
#     for trait in traits
#         # trait = traits[1]
#         @show trait
#         genomes, phenomes = begin
#             y = df[!, trait]
#             idx = findall(.!ismissing.(y) .&& .!isinf.(y) .&& .!isnan.(y))
#             P_0 = Phenomes(n=length(idx), t=1)
#             P_0.entries = df.entries[idx]
#             P_0.populations .= "bulked"
#             P_0.traits[1] = trait
#             P_0.phenotypes[:, 1] = y[idx]
#             merge(G_0, P_0)
#         end
#         fit = GenomicBreedingModels.bayesa(genomes=genomes, phenomes=phenomes)
#         fits[trait] = fit
#         # fit.b_hat
#         # fit.metrics["h²"]
#         # h²_realised[trait]
#     end
#     JLD2.save("fits.jld2", fits)
#     fits
# end

# Genomic prediction models
function gpmodels(
    genomes::Genomes,
    phenomes::Phenomes;
    save::Bool=true,
)::Dict{String,Fit}
    fits = Dict()
    for trait in phenomes.traits
        # trait = phenomes.traits[1]
        @show trait
        G, P = begin
            idx_traits = findall(phenomes.traits .== trait)
            y = phenomes.phenotypes[:, idx_traits][:, 1]
            idx_entries = findall(.!ismissing.(y) .&& .!isinf.(y) .&& .!isnan.(y))
            P = slice(phenomes, idx_entries=idx_entries, idx_traits=idx_traits)
            merge(genomes, P)
        end
        fit = GenomicBreedingModels.bayesa(genomes=G, phenomes=P)
        fits[trait] = fit
    end
    if save
        JLD2.save("fits.jld2", fits)
    end
    fits
end


# GEBVs
function gebvs(
    fit::Fit,
    genomes::Genomes;
    idx_entries::Union{Nothing,Vector{Int64}}=nothing,
)::Phenomes
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        Int64.(idx_entries)
    end
    phenomes = Phenomes(n=length(idx_entries), t=1)
    phenomes.entries = genomes.entries
    phenomes.populations = genomes.populations
    phenomes.phenotypes[:, 1] = GenomicBreedingModels.predict(fit=fit, genomes=genomes, idx_entries=idx_entries)
    phenomes
end


# Simulate mating
function panmixis(
    genomes::Genomes;
    idx_entries::Union{Nothing,Vector{Int64}}=nothing,
    population_name::String="",
    n::Int64=1_000,
    seed::Int64=42,
    verbose::Bool=false
)::Genomes
    # genomes=GenomicBreedingCore.simulategenomes(n=119, l=10_000); idx_entries=nothing; n=1_000; seed=42
    Random.seed!(seed)
    idx_entries = if isnothing(idx_entries)
        collect(1:length(genomes.entries))
    else
        Int64.(idx_entries)
    end
    idx_0 = sample(idx_entries, n)
    idx_1 = sample(idx_entries, n)
    progenies = Genomes(n=n, p=length(genomes.loci_alleles))
    progenies.entries = ["entry-$(idx_0[i])x$(idx_1[i])-$(rand(1:n))" for i in 1:n]
    progenies.populations .= population_name
    progenies.loci_alleles = genomes.loci_alleles
    if verbose
        pb = ProgressMeter.Progress(n, "Simulating random mating")
    end
    for (i, (j, k)) in enumerate(zip(idx_0, idx_1))
        # i = 1; j = idx_0[i]; k = idx_1[i]
        progenies.allele_frequencies[i, :] = (genomes.allele_frequencies[j, :] .+ genomes.allele_frequencies[k, :]) ./ 2
        if verbose
            ProgressMeter.next!()
        end
    end
    if verbose
        ProgressMeter.finish!()
    end
    progenies
end


# seed = 42
# n_cycles = 3
# population_size = 1_000
# selection_intensity = 0.1
# data = Dict()
# for trait in traits
#     # trait = traits[1]
#     @show trait
#     fit = fits[trait]
#     genomes = clone(G_0)
#     phenomes = gebvs(; fit=fit, genomes=genomes)
#     data[trait] = Dict(
#         "phenomes" => [phenomes], 
#         "D" => [Normal(mean(phenomes.phenotypes[:,1]), std(phenomes.phenotypes[:,1]))]
#     )
#     for t in 1:n_cycles
#         # t = 1
#         idx_entries = begin
#             idx = sortperm(phenomes.phenotypes[:, 1], rev=true)
#             n = length(idx)
#             s = Int64(ceil(selection_intensity * n))
#             idx[1:s]
#         end
#         selected_genomes = slice(genomes, idx_entries=idx_entries)
#         @assert checkdims(selected_genomes)
#         genomes = panmixis(selected_genomes, n=population_size, population_name="cycle_$t", seed=seed)
#         @assert checkdims(genomes)
#         phenomes = gebvs(; fit=fit, genomes=genomes)
#         push!(data[trait]["phenomes"], phenomes)
#         push!(data[trait]["D"], Normal(mean(phenomes.phenotypes[:,1]), std(phenomes.phenotypes[:,1])))
#     end
# end
# JLD2.save("data.jld2", data)

function multigenphenomes(
    genomes::Genomes,
    phenomes::Phenomes;
    fits::Dict{String,Fit},
    n_cycles::Int64=3,
    population_size::Int64=1_000,
    selection_intensity::Float64=0.1,
    save::Bool=true,
    seed::Int64=42,
)::Dict{String,Union{Phenomes,Normal}}
    data = Dict()
    for trait in phenomes.traits
        # trait = traits[1]
        @show trait
        fit = fits[trait]
        G = clone(genomes)
        P = gebvs(; fit=fit, genomes=G)
        data[trait] = Dict(
            "phenomes" => [P],
            "D" => [Normal(mean(P.phenotypes[:, 1]), std(P.phenotypes[:, 1]))]
        )
        for t in 1:n_cycles
            # t = 1
            idx_entries = begin
                idx = sortperm(P.phenotypes[:, 1], rev=true)
                n = length(idx)
                s = Int64(ceil(selection_intensity * n))
                idx[1:s]
            end
            selected_G = slice(G, idx_entries=idx_entries)
            @assert checkdims(selected_G)
            G = panmixis(selected_G, n=population_size, population_name="cycle_$t", seed=seed)
            @assert checkdims(G)
            P = gebvs(fit, G)
            push!(data[trait]["phenomes"], P)
            push!(data[trait]["D"], Normal(mean(P.phenotypes[:, 1]), std(P.phenotypes[:, 1])))
        end
    end
    if save
        JLD2.save("data.jld2", data)
    end
    data
end


# Unicodeplots of genetic gain distributions
for trait in traits
    # trait = traits[1]
    @show trait
    y, d = begin
        y = sort(rand(data[trait]["D"][1], population_size))
        d = pdf(data[trait]["D"][1], y)
        d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
        (y, d)
    end
    xlim = ((minimum(y)), maximum(vcat(y, data[trait]["phenomes"][end].phenotypes[:, 1])))
    p = UnicodePlots.lineplot(y, d, title="$trait (h²=$(round(fits[trait].metrics["h²"], digits=2)))", width=100, xlim=xlim)
    for i in 2:length(data[trait]["phenomes"])
        y, d = begin
            y = sort(rand(data[trait]["D"][i], population_size))
            d = pdf(data[trait]["D"][i], y)
            d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
            (y, d)
        end
        UnicodePlots.lineplot!(p, y, d)
    end
    display(p)
end
# Plot genomically informed genetic gain distributions and save as SVG
fig = CairoMakie.Figure(size=(500, 1_500))
N = length(traits)
@assert N == 14
p = 7;
q = 2;
axs = Dict()
for i in 0:(p-1)
    for j in 0:(q-1)
        # i = 0; j = 0        
        idx = (i * q) + j
        trait = traits[idx+1]
        y, d = begin
            y = sort(rand(data[trait]["D"][1], population_size))
            d = pdf(data[trait]["D"][1], y)
            d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
            (y, d)
        end
        xlim = ((minimum(y)), maximum(vcat(y, data[trait]["phenomes"][end].phenotypes[:, 1])))
        title = string(split(trait, "|")[2], "\n", replace(replace(replace(split(trait, "|")[3], "NUE" => "highN"), "WUE" => "drought"), "Control" => "lowN"))
        axs[trait] = CairoMakie.Axis(
            fig[i, j],
            title="$title\n(h²=$(round(fits[trait].metrics["h²"],digits=2)))",
            limits=(xlim, nothing)
        )
        CairoMakie.lines!(axs[trait], y, d)
    end
end
# Plot
for trait in traits
    # trait = traits[1]
    for i in 2:length(data[trait]["phenomes"])
        y, d = begin
            y = sort(rand(data[trait]["D"][i], population_size))
            d = pdf(data[trait]["D"][i], y)
            d = (d .- minimum(d)) ./ (maximum(d) - minimum(d))
            (y, d)
        end
        CairoMakie.lines!(axs[trait], y, d)
    end
end
CairoMakie.save("sample_selection_simulations_less_naive.svg", fig)

