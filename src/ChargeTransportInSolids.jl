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
export tiny_penalty_value

##################################################################

include("ct_distributions.jl")

export Boltzmann, Blakemore, FermiDiracMinusOne, FermiDiracOneHalfBednarczyk, FermiDiracOneHalfTeSCA, FermiDiracZero

##################################################################

include("ct_datatypes.jl")

export grid_dimension, OneD_grid, TwoD_grid, ThreeD_gride

export boundary_model, ohmic_contact, schottky_contact
export interface_model, interface_model_none, interface_model_surface_recombination, interface_model_ion_charge, interface_model_tangential_flux, interface_model_surface_recombination_and_tangential_flux

export interface_model_discont_qF

export model_type, model_transient, model_stationary

export abstract_SRH_model, SRH_model, SRH_model_off, SRH_model_without_traps, SRH_model_without_traps_stationary, SRH_model_with_traps,  SRH_model_traps_transient

export SRH_2species_present_trap_dens
##################################################################

include("ct_physics.jl")

export ScharfetterGummel, excessChemicalPotential, diffusionEnhanced, generalized_SG 
export inEquilibrium, outOfEquilibrium
export ScharfetterGummel_Graded, excessChemicalPotential_Graded
export generation_model, generation_none, generation_beer_lambert, generation_uniform
export breaction!, bstorage!, reaction!, storage!, flux!

##################################################################

include("ct_system.jl")

export ChargeTransportParams, ChargeTransportParamsNodal, ChargeTransportData, ChargeTransportSystem
export ChargeTransportBulkRecombination, set_bulk_recombination, ChargeTransportIonicChargeCarriers, enable_ion_vacancies
export equilibrium_solve!
export set_ohmic_contact!, set_schottky_contact!
export compute_densities!, compute_energies!, electroNeutralSolution!, print_jacobi
export show_params, trap_density!
export linearScanProtocol, scan_protocol_type, set_time_mesh, get_current_val
export chargeDensity
export enable_traps!

##################################################################

include("ct_plotting.jl")

export plot_densities, plot_energies, plot_doping, plot_electroNeutralSolutionBoltzmann, plot_solution, plot_IV


end # module
