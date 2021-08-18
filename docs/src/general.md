
ChargeTransportInSolids.jl -- A drift diffusion solver 
================================
This package is a prototype for solving drift-diffusion equations for the simulation of charge transport in solids that centers physics preserving discretization schemes.

Concerning the spatial discretization we rely on  a finite volume method implemented within the solver [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl). Currently, we use for the time discretization an implicit Euler method.

!!! note

    The package exclusively assumes charge transport portrayed by drift-diffusion equations (Poisson equation + continuity equation(s)).




Installation and First Steps
================================
The installation can be easily done via the Julia REPL by the following commands

```julia
julia> using Pkg
julia> Pkg.add("ChargeTransportInSolids")

## or within the package environment
julia> ]
(@v1.6) pkg> add ChargeTransportInSolids
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
We recommend on having a look either on the Pluto notebook or on the example files. The examples can be loaded and run by 

```julia
julia> include("Example103_PSC.jl")
Main.Example103_PSC
julia> Example103_PSC.main()
```
Due to the encapsulation of the examples into modules, you can load as many examples as you like. If you want to modify one of the examples, consider using [Revise.jl](https://github.com/timholy/Revise.jl) and `includet`.

A guide on how to understand the main parts of the example problems and, thus, how to properly build your own class of problems is explained here:

```@contents
Pages = [
    "examples/example-GaAs.md",
    "examples/example-PSC.md",
    ]
Depth = 2
```
