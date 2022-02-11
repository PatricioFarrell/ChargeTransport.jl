module ChargeTransport

using VoronoiFVM
using ExtendableGrids
using Printf
using DocStringExtensions
using SparseArrays
using Roots


include("ct_constants.jl")

export kB, Planck_constant, mₑ, q, ε0


include("ct_units.jl")

export K, J, A, V, m, s, C, kg
export cm, mm, μm, nm, ms, μs, ns, ps, eV
export tiny_penalty_value

##################################################################

include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk
export FermiDiracOneHalfTeSCA, FermiDiracZero

##################################################################

include("ct_datatypes.jl")

export BoundaryModel, OhmicContact, SchottkyContact
export InterfaceModel, InterfaceModelNone, InterfaceModelSurfaceReco
export InterfaceModelIonCharge, InterfaceModelTangentialFlux
export InterfaceModelSurfaceRecoAndTangentialFlux

export InterfaceModelDiscontqF

export ModelType, Transient, Stationary

export SRHOff, SRHWithoutTraps
export SRHWithoutTrapsStationary, SRHWithTraps, SRHTrapsTransient

export InEquilibrium, OutOfEquilibrium
##################################################################

include("ct_physics.jl")

export ScharfetterGummel, ExcessChemicalPotential, DiffusionEnhanced, GeneralizedSG
export ScharfetterGummelGraded, ExcessChemicalPotentialGraded
export GenerationModel, GenerationNone, GenerationBeerLambert, GenerationUniform
export breaction!, bstorage!, reaction!, storage!, flux!

##################################################################

include("ct_system.jl")

export Params, ParamsNodal, Data, System
export BulkRecombination, set_bulk_recombination, IonicChargeCarriers, enable_ionic_carriers
export equilibrium_solve!
export set_contact!
export compute_densities!, compute_energies!, electroNeutralSolution!, print_jacobi
export show_params, trap_density!
export LinearScanProtocol, ScanProtocolType, set_time_mesh, get_current_val
export charge_density
export enable_traps!

##################################################################

include("ct_plotting.jl")

export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann
export plot_solution, plot_IV


end # module
