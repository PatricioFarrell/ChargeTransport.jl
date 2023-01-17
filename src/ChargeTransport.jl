module ChargeTransport

using VoronoiFVM
using ExtendableGrids
using Printf
using DocStringExtensions
using SparseArrays
using Interpolations
using Roots


include("ct_constants.jl")

export kB, Planck_constant, mₑ, q, ε0
##################################################################

include("ct_units.jl")

export K, J, A, V, m, s, C, kg, Hz, kHz, W, kW
export cm, mm, μm, nm, ms, μs, ns, ps, eV
export tiny_penalty_value
##################################################################

include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk
export FermiDiracOneHalfTeSCA, FermiDiracZero
##################################################################

include("ct_datatypes.jl")

export StandardFuncSet
export QType
export QFModelType, DiscontQF, ContQF

export OuterBoundaryModelType, OuterBoundaryModelType, InterfaceModelType
export OhmicContact, SchottkyContact, SchottkyBarrierLowering
export InterfaceNone, InterfaceRecombination

export ModelType, Transient, Stationary

export FluxApproximationType
export ScharfetterGummel, ExcessChemicalPotential, DiffusionEnhanced, DiffusionEnhancedModifiedDrift, GeneralizedSG
export ScharfetterGummelGraded, ExcessChemicalPotentialGraded

export InEquilibrium, OutOfEquilibrium

export SRHModelType, SRHWithoutTrapsType, SRHWithTrapsType
export SRHOff, SRHWithoutTrapsStationary, SRHTrapsTransient

export GenerationModelType
export GenerationNone, GenerationBeerLambert, GenerationUniform
export BarrierLoweringType
export BarrierLoweringOn, BarrierLoweringOff
##################################################################

include("ct_physics.jl")

export breaction!, bstorage!, reaction!, storage!, flux!
export zeroVoltage
##################################################################

include("ct_system.jl")

export Params, ParamsNodal, Data, System
export BulkRecombination, set_bulk_recombination

export enable_ionic_carrier!, enable_trap_carrier!

export equilibrium_solve!
export set_contact!
export compute_densities!, compute_energies!
export compute_open_circuit_voltage
export electroNeutralSolution!, print_jacobi
export show_params, trap_density!
export get_current_val, charge_density

##################################################################

include("ct_plotting.jl")

export set_plotting_labels
export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann
export plot_solution, plot_IV


end # module