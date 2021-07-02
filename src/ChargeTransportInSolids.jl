module ChargeTransportInSolids

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


include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA

include("ct_physics.jl")

export boundary_model, ohmic_contact, schottky_contact
export interface_model, interface_model_none, interface_model_surface_recombination, interface_model_ion_charge, interface_model_ion_charge_left, interface_model_ion_charge_right
export model_type, model_transient, model_stationary
export bulk_recombination_model, bulk_recombination_none, bulk_recombination_trap_assisted, bulk_recombination_radiative, bulk_recombination_full
export ScharfetterGummel, excessChemicalPotential, diffusionEnhanced, generalized_SG 
export inEquilibrium, outOfEquilibrium
export ScharfetterGummel_Graded, excessChemicalPotential_Graded
export breaction!, bstorage!, reaction!, storage!, flux!

include("ct_system.jl")

export ChargeTransportParams, ChargeTransportParamsNodal, ChargeTransportData, ChargeTransportSystem
export ct_enable_species!, ct_enable_boundary_species!, ct_unknowns, ct_solve!
export set_ohmic_contact!, set_indices!, compute_densities!, compute_energies!, electroNeutralSolution!, print_jacobi
export trap_density!

include("ct_plotting.jl")

export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann, plot_solution, plot_IV



end # module
