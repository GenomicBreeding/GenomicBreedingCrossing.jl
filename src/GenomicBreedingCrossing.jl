module GenomicBreedingCrossing

using GenomicBreedingCore, GenomicBreedingModels
using Random, StatsBase, DataFrames, Distributions, UnicodePlots, MultivariateStats, JLD2
using CairoMakie, ColorSchemes

include("sim_naive.jl")
include("sim_gebv.jl")

end
