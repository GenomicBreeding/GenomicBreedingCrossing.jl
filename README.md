# GenomicBreedingCrossing

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GenomicBreeding.github.io/GenomicBreedingCrossing.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GenomicBreeding.github.io/GenomicBreedingCrossing.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreedingCrossing.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GenomicBreeding/GenomicBreedingCrossing.jl/actions/workflows/CI.yml?query=branch%3Amain)

Cross/mating designs and breeding simulations. It provides tools for defining breeding populations, performing crosses, fitting genomic prediction models, and simulating progeny with and without genome-estimated breeding values.

## Dev stuff:

### REPL prelude

```shell
julia --project=. --threads=2,1 --load test/interactive_prelude.jl
```

### Format and test

```shell
time julia --project=. --threads=2 -e "using Pkg; Pkg.update()"
time julia --project=. --threads=2  test/cli_tester.jl
```
