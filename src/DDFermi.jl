module DDFermi

using VoronoiFVM
using PyPlot; PyPlot.pygui(true)
using Printf
using DocStringExtensions

include("dd_constants.jl")
include("dd_units.jl")
include("dd_distributions.jl")
include("dd_system.jl")

# export K, J, A, V, m, s, C, kg
# export cm, mm, μm, nm, eV
# export kB, Planck_constant, mₑ, q, ε0
# export Boltzmann, Blakemore, FermiDirac

end # module
