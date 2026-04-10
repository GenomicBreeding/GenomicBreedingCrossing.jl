module GenomicBreedingCrossing

using GenomicBreedingCore, GenomicBreedingModels
using Random,
    UUIDs,
    Dates,
    StatsBase,
    DataFrames,
    Distributions,
    UnicodePlots,
    MultivariateStats,
    LinearAlgebra,
    JLD2,
    CSV,
    ProgressMeter
using CairoMakie, ColorSchemes

include("breeding_populations.jl")
export BreedingPopulations, check_populations, plot_simple, plot_and_save

include("sim_naive.jl")
export simulate_naive

include("crossing.jl")
export panmixis, paired_cross

include("sim_gebv.jl")
export fit_gp_models, predict_gebvs, simulate_gebv

include("optim.jl")
# export ...

end
