ChargeTransportInSolids.jl -- A drift-diffusion solver 
================================

[![Build status](https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/workflows/linux-macos-windows/badge.svg)](https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PatricioFarrell.github.io/ChargeTransportInSolids.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PatricioFarrell.github.io/ChargeTransportInSolids.jl/dev)


This package is a prototype for solving drift-diffusion equations for the simulation of charge transport in solids that centers physics preserving discretization schemes.

Concerning the spatial discretization we rely on a Voronoi finite volume method implemented within the solver [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl). Currently, we use for the time discretization an implicit Euler method.

!!! note

    This package exclusively assumes charge transport portrayed by drift-diffusion equations (Poisson equation + continuity equation(s)).


The simulations in following papers are based on ChargeTransportInSolids.jl

[1.] D. Abdel, P. Farrell and J. Fuhrmann.[ Assessing the quality of the excess chemical potential flux scheme for degenerate semiconductor device simulation.](https://link.springer.com/article/10.1007/s11082-021-02803-4) Optical and Quantum Electronics 53 (163) (2021).

[2.] D. Abdel, P. Vágner, J. Fuhrmann and P. Farrell. [Modelling charge transport in perovskite solar cells: Potential-based and limiting ion depletion.](https://www.sciencedirect.com/science/article/abs/pii/S0013468621009865) Electrochimica Acta 390 (2021).

Installation and First Steps
================================
The installation can be easily done via the Julia REPL by the following commands

```julia
julia> using Pkg
julia> Pkg.add("ChargeTransportInSolids")

## or within the package environment
julia> ]
(@v1.6.1) pkg> add ChargeTransportInSolids
```

For the construction of an example the following packages need to be loaded

```julia
# concrete application dependent numerical tools
julia> using ChargeTransportInSolids
# nonlinear partial differential equations solver (based on a Voronoi finite volume method)
julia> using VoronoiFVM
# package for storing grid information
julia> using ExtendableGrids
```
We recommend on having a look on the example files. A guide on how to understand the main parts of the example problems and, thus, how to properly build your own class of problems is explained here:

```@contents
Pages = [
    "examples/GaAs.md",
    "examples/PSC.md",
    ]
Depth = 2
```

Further, the examples can be loaded and run by 

```julia
julia> include("Example103_PSC.jl")
Main.Example103_PSC
julia> Example103_PSC.main()

## or if you are interested in visualization of solutions ## (currently, predefined functions only tested with PyPlot)
julia> include("Example103_PSC.jl")
Main.Example103_PSC
julia> Example103_PSC.main(plotting = true)
```
Due to the encapsulation of the examples into modules, you can load as many examples as you like. If you want to modify one of the examples, consider using [Revise.jl](https://github.com/timholy/Revise.jl) and `includet`.

