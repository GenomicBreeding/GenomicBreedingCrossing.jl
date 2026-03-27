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
    JLD2,
    CSV,
    ProgressMeter
using CairoMakie, ColorSchemes

include("breeding_populations.jl")
export BreedingPopulations, check_populations, plot_simple, plot_and_save

include("sim_naive.jl")
export simulate_naive

include("sim_gebv.jl")
export fit_gp_models, simulate_gebv

end
