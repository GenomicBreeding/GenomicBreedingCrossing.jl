"""
    panmixis(genomes::Genomes; idx_entries::Union{Nothing,Vector{Int64}} = nothing, 
             population_name::String = "", n::Int64 = 1_000, seed::Int64 = 42, 
             ϵ::Float64 = 0.001, verbose::Bool = false)::Genomes

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
- `simple_mating_model::Bool`: Use the simplistic allele frequency mean crossing model, otherwise
   simulate linkage disequillibrium for use in sampling from a multivariate normal distribution. Default: `false`
- `seed::Int64`: Random seed for reproducibility. Default: `42`
- `ϵ::Float64`: Value used for the diagonal inflation of the LD/covariance matrix, i.e. when `simple_mating_model=false`. Default: `0.001`
- `verbose::Bool`: If `true`, display progress bar during simulation. Default: `false`

# Returns
- `Genomes`: A new `Genomes` object containing the simulated progeny with:
  - Allele frequencies computed as the mean of two randomly sampled parents
  - Unique UUIDs for each progeny entry
  - Same loci structure as the source genomes
  - Assigned population name

# Examples
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
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
    simple_mating_model::Bool = false,
    ϵ::Float64 = 0.001,
    seed::Int64 = 42,
    verbose::Bool = false,
)::Genomes
    # genomes = GenomicBreedingCore.simulategenomes(n=119, l=10_000); idx_entries=nothing; population_name::String = ""; n=1_000; simple_mating_model=false; seed=42; ϵ = 0.001; verbose = true
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
    if simple_mating_model
        # Very simplistic approach assuming genotypes are pools and LD is absent!
        pb = verbose ? ProgressMeter.Progress(n, "Simulating random mating") : nothing
        for (i, (j, k)) in enumerate(zip(idx_0, idx_1))
            # i = 1; j = idx_0[i]; k = idx_1[i]
            progenies.allele_frequencies[i, :] =
                (genomes.allele_frequencies[j, :] .+ genomes.allele_frequencies[k, :]) ./ 2
            verbose ? ProgressMeter.next!(pb) : nothing
        end
        verbose ? ProgressMeter.finish!(pb) : nothing
    else
        verbose ? println("Simulating random mating") : nothing
        verbose ? println("Estimating LD") : nothing
        # Estimate linkage disequillibrium 
        Σ = cov(genomes.allele_frequencies)
        verbose ? println("Making sure the LD matrix is positive semi-definite") : nothing
        # Set the variance/covariance of fixed loci to zero
        Σ[isinf.(Σ)] .= 0.0
        Σ[isnan.(Σ)] .= 0.0
        counter = 0
        while !isposdef(Σ) && (counter < 10)
            Σ[diagind(Σ)] .+= ϵ
            counter += 1
        end
        if !isposdef(Σ)
            throw(
                ErrorException(
                    "Uggghhhh... Cannot estimate the LD matrix! Send me an email to implement a proper fix for this. Sorry :-(",
                ),
            )
        end
        verbose ? println("Extracting allele frequency means") : nothing
        μ = Float64.(mean(genomes.allele_frequencies, dims = 1)[1, :])
        verbose ?
        println(
            "Defining the multivariate distribution using the means of allele frequencies and LD matrix",
        ) : nothing
        D = Distributions.MvNormal(μ, Σ)
        verbose ? println("Sampling new genotypes from the multivariate distribution") :
        nothing
        G = Matrix(rand(rng, D, n)')
        verbose ? println("Mapping the sampled values into the zero to one range") : nothing
        G[G .< 0.0] .= 0.0
        G[G .> 1.0] .= 1.0
        progenies.allele_frequencies = G
        # if verbose
        #     # Commenting these out as correlation computations can be very expensive
        #     display(UnicodePlots.heatmap(first(values(dist)), title="Parents' LD matrix"))
        #     display(UnicodePlots.heatmap(cor(G), title="Progenies' LD matrix"))
        # end
    end
    @assert checkdims(progenies) # TODO: emit a better more specific error message here
    progenies
end

"""
    paired_cross(
        genomes_1::Genomes,
        genomes_2::Genomes; Σ::Union{Nothing,
        Matrix{Float64}} = nothing,
        population_name::String = "",
        n::Int64 = 1_000,
        discrete_loci_alleles::Bool = false,
        ϵ::Float64 = 0.001,
        seed::Int64 = 42,
        verbose::Bool = false
    )::Genomes

Simulate paired crossing between two parent genomes to generate progeny.

This function performs crossing between two specific parent genomes, creating progeny with allele frequencies sampled from a multivariate normal distribution based on the parents' allele frequencies and linkage disequilibrium.

# Arguments
- `genomes_1::Genomes`: First parent genome (must contain exactly one entry)
- `genomes_2::Genomes`: Second parent genome (must contain exactly one entry)
- `Σ::Union{Nothing, Matrix{Float64}}`: Covariance matrix for the multivariate distribution. If `nothing`, uses identity matrix. Default: `nothing`
- `population_name::String`: Name to assign to the progeny population. Default: `""`
- `n::Int64`: Number of progeny to generate. Default: `1_000`
- `discrete_loci_alleles::Bool`: If `true`, discretize allele frequencies to match parent values. Default: `false`
- `ϵ::Float64`: Value used for the diagonal inflation of the LD/covariance matrix, i.e. when `simple_mating_model=false`. Default: `0.001`
- `seed::Int64`: Random seed for reproducibility. Default: `42`
- `verbose::Bool`: If `true`, print progress information. Default: `false`

# Returns
- `Genomes`: A new `Genomes` object containing the simulated progeny

# Examples: TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingModels, GenomicBreedingCrossing, StatsBase, DataFrames, Combinatorics)
julia> genomes = GenomicBreedingCore.simulategenomes(n=119, l=10_000, verbose=false);

julia> genomes_1 = slice(genomes, idx_entries=[1]);

julia> genomes_2 = slice(genomes, idx_entries=[2]);

julia> progenies = paired_cross(genomes_1, genomes_2, Σ=cov(genomes.allele_frequencies), n=123);

julia> size(progenies.allele_frequencies)
(123, 10000)

julia> isnothing(progenies.allele_frequencies_homologous_chroms)
true

julia> genomes.allele_frequencies_homologous_chroms = 1.00 .- genomes.allele_frequencies;

julia> genomes_1 = slice(genomes, idx_entries=[1]);

julia> genomes_2 = slice(genomes, idx_entries=[2]);

julia> progenies = paired_cross(genomes_1, genomes_2, Σ=cov(genomes.allele_frequencies), n=123);

julia> size(progenies.allele_frequencies)
(123, 10000)

julia> isnothing(progenies.allele_frequencies_homologous_chroms)
false
```
"""
function paired_cross(
    genomes_1::Genomes,
    genomes_2::Genomes;
    Σ::Union{Nothing,Matrix{Float64}} = nothing,
    population_name::String = "",
    n::Int64 = 1_000,
    seed::Int64 = 42,
    ϵ::Float64 = 0.001,
    verbose::Bool = false,
)::Genomes
    # genomes = GenomicBreedingCore.simulategenomes(n=119, l=10_000); population_name::String = ""; n=1_000; seed=42; ϵ::Float64 = 0.001; verbose = true
    # genomes.allele_frequencies_homologous_chroms = 1.00 .- genomes.allele_frequencies
    # genomes_1 = slice(genomes, idx_entries=[1])
    # genomes_2 = slice(genomes, idx_entries=[2])
    # Σ = cov(genomes.allele_frequencies)
    # Σ[isinf.(Σ)] .= 0.0
    # Σ[isnan.(Σ)] .= 0.0
    # counter = 0
    # while !isposdef(Σ) && (counter < 10)
    #     Σ[diagind(Σ)] .+= ϵ
    #     counter += 1
    # end
    ##############################################################
    if (length(genomes_1.entries) == length(genomes_2.entries) == 1) == false
        throw(ErrorException("We expect 1 entry in each genome for paired crossing"))
    end
    if (length(genomes_1.loci_alleles) == length(genomes_2.loci_alleles)) == false
        throw(
            ErrorException(
                "Incompatible number of loci-alleles: $(length(genomes_1.loci_alleles)) and $(length(genomes_2.loci_alleles))",
            ),
        )
    end
    Random.seed!(seed)
    rng = Xoshiro(seed)
    _, p = size(genomes_1.allele_frequencies)
    Σ::Matrix{Float64} = if isnothing(Σ)
        I = zeros(p, p)
        I[diagind(i)] .= 1.00
        I
    else
        Σ
    end
    (G_1, G_2) =
        if isnothing(genomes_1.allele_frequencies_homologous_chroms) ||
           isnothing(genomes_2.allele_frequencies_homologous_chroms)
            ########################################
            ### Does not use phasing information ###
            ########################################
            μ =
                (genomes_1.allele_frequencies[1, :] .+ genomes_2.allele_frequencies[1, :]) ./
                2
            verbose ? println("Making sure the LD matrix is positive semi-definite") :
            nothing
            # Set the variance/covariance of fixed loci to zero
            Σ[isinf.(Σ)] .= 0.0
            Σ[isnan.(Σ)] .= 0.0
            counter = 0
            while !isposdef(Σ) && (counter < 10)
                Σ[diagind(Σ)] .+= ϵ
                counter += 1
            end
            if !isposdef(Σ)
                throw(
                    ErrorException(
                        "Uggghhhh... Cannot estimate the LD matrix! Send me an email to implement a proper fix for this. Sorry :-(",
                    ),
                )
            end
            D = Distributions.MvNormal(μ, Σ)
            G = Matrix(rand(rng, D, n)')
            δ_1 = abs.(G .- genomes_1.allele_frequencies)
            δ_2 = abs.(G .- genomes_2.allele_frequencies)
            idx_1 = δ_1 .<= δ_2
            idx_2 = δ_1 .> δ_2
            for j = 1:p
                # j = 1
                G[idx_1[:, j], j] .= genomes_1.allele_frequencies[1, j]
                G[idx_2[:, j], j] .= genomes_2.allele_frequencies[1, j]
            end
            (G, nothing)
        else
            ################################
            ### Uses phasing information ###
            ################################
            # Instantiate the offspring genomes at the first locus-allele position
            G_1 = zeros(n, p)
            G_2 = 1.00 .- G_1
            G_1[:, 1] = sample(
                rng,
                [
                    genomes_1.allele_frequencies[1, 1],
                    genomes_1.allele_frequencies_homologous_chroms[1, 1],
                ],
                n,
                replace = true,
            )
            G_2[:, 1] = sample(
                rng,
                [
                    genomes_2.allele_frequencies[1, 1],
                    genomes_2.allele_frequencies_homologous_chroms[1, 1],
                ],
                n,
                replace = true,
            )
            # Define the likelihood linkage per pair of loci
            for j = 2:p
                # j = 2
                v_prev = Σ[(j-1), (j-1)]
                v_curr = Σ[j, j]
                corr = Σ[(j-1), j] / (sqrt(v_prev) * sqrt(v_curr))
                # From genomes_1, i.e. parent 1
                let
                    bool_linkage = rand(rng, n) .<= corr
                    bool_homologue_A = G_1[:, j-1] .== genomes_1.allele_frequencies[1, j-1]
                    bool_homologue_B =
                        G_1[:, j-1] .==
                        genomes_1.allele_frequencies_homologous_chroms[1, j-1]
                    G_1[bool_homologue_A.&&bool_linkage, j] .=
                        genomes_1.allele_frequencies[1, j]
                    G_1[bool_homologue_B.&&bool_linkage, j] .=
                        genomes_1.allele_frequencies_homologous_chroms[1, j]
                    G_1[bool_homologue_A.&&.!bool_linkage, j] .=
                        genomes_1.allele_frequencies_homologous_chroms[1, j]
                    G_1[bool_homologue_B.&&.!bool_linkage, j] .= genomes_1.allele_frequencies[1, j]
                    @assert n ==
                            sum(bool_homologue_A .&& bool_linkage) +
                            sum(bool_homologue_B .&& bool_linkage) +
                            sum(bool_homologue_A .&& .!bool_linkage) +
                            sum(bool_homologue_B .&& .!bool_linkage)
                end
                # From genomes_2, i.e. parent 2
                let
                    bool_linkage = rand(rng, n) .<= corr
                    bool_homologue_A = G_2[:, j-1] .== genomes_2.allele_frequencies[1, j-1]
                    bool_homologue_B =
                        G_2[:, j-1] .==
                        genomes_2.allele_frequencies_homologous_chroms[1, j-1]
                    G_2[bool_homologue_A.&&bool_linkage, j] .=
                        genomes_2.allele_frequencies[1, j]
                    G_2[bool_homologue_B.&&bool_linkage, j] .=
                        genomes_2.allele_frequencies_homologous_chroms[1, j]
                    G_2[bool_homologue_A.&&.!bool_linkage, j] .=
                        genomes_2.allele_frequencies_homologous_chroms[1, j]
                    G_2[bool_homologue_B.&&.!bool_linkage, j] .= genomes_2.allele_frequencies[1, j]
                    @assert n ==
                            sum(bool_homologue_A .&& bool_linkage) +
                            sum(bool_homologue_B .&& bool_linkage) +
                            sum(bool_homologue_A .&& .!bool_linkage) +
                            sum(bool_homologue_B .&& .!bool_linkage)
                end
            end
            (G_1, G_2)
        end
    progenies = Genomes(n = n, p = length(genomes_1.loci_alleles))
    ids = [uuid4(rng) for i = 1:n]
    progenies.entries = string.("entry-", ids)
    progenies.populations .= population_name
    progenies.loci_alleles = genomes_1.loci_alleles
    progenies.allele_frequencies = G_1
    progenies.allele_frequencies_homologous_chroms = G_2
    # if verbose
    #     # Commenting these out as correlation computations can be very expensive
    #     display(UnicodePlots.heatmap(Σ, title="Parents' LD matrix"))
    #     display(UnicodePlots.heatmap(cov(G_1), title="Progenies' LD matrix"))
    # end
    @assert checkdims(progenies) # TODO: emit a better more specific error message here
    progenies
end
