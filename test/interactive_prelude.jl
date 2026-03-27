using Pkg
Pkg.activate(".")
using GenomicBreedingCrossing
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
