Running the examples
====================
JUST COPY AND PASTED FROM VORONOIFVM!!!!

The examples have been designed with the following issues in mind:
- they run from the Julia REPL
- each example is a Julia module named similar to the basename of the example file.
- an example can be uses as the starting point for a project 
- the examples at the same time comprise the test suite for VoronoiFVM.

Plotting is performed using the GridVisualize.jl package which interfaces PyPlot.jl,
Plots.jl, Makie.jl.

In order to run `ExampleXXX`, peform the following steps:

- Download the example file (e.g. via the source code link at the top)
- Call Julia with  an Julia environment which contains `VoronoiFVM.jl`, `ExtendableGrids.jl`, `GridVisualize.jl` 
  and e.g. `PyPlot.jl`
- `include("ExampleXXX.jl")`
- Run the example via `ExampleXXX.main(Plotter=PyPlot)`

Due to the encapsulation into modules, you can load as many examples as you like.

If you want to modifiy the example, consider using `Revise.jl` and `includet`.  
