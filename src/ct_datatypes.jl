##########################################################
##########################################################
"""

$(TYPEDEF)
Abstract type for boundary model. Subtypes are ohmic_contact, schottky_contact
and interface_model.

"""
abstract type boundary_model   end

############    outer boundary conditions     ############
"""
Abstract type for ohmic contacts as outer boundary model.

"""
abstract type ohmic_contact <: boundary_model  end


"""
Abstract type for schottky contacts as boundary model.

"""
abstract type schottky_contact <: boundary_model end

############    inner boundary conditions     ############
"""
$(TYPEDEF)
Abstract type for interface model which
is part of boundary model with several subtypes.

"""
abstract type interface_model <: boundary_model end

abstract type interface_model_tangential_flux                           <: interface_model end

abstract type interface_model_surface_recombination_and_tangential_flux <: interface_model end

"""
$(TYPEDEF)
Abstract type for no interface model.

"""
abstract type interface_model_none <: interface_model end

#
#$(TYPEDEF)
#Abstract type for an interface model where discontinuous
#quasi Fermi potentials are needed.
#
abstract type interface_model_discont_qF <: interface_model end


"""
$(TYPEDEF)
Abstract type for surface recombination mechanisms.

"""
abstract type interface_model_surface_recombination <: interface_model end


# """
# $(TYPEDEF)
# Abstract type for present ion charges at interfaces.

# """
abstract type interface_model_ion_charge <: interface_model end

##########################################################
##########################################################

"""
$(TYPEDEF)
model_SRH_without_traps as parent of different subtypes.

"""
abstract type model_SRH_without_traps end

abstract type model_SRH_stationary     <: model_SRH_without_traps end
abstract type model_SRH_off            <: model_SRH_without_traps end

"""
$(TYPEDEF)
model_SRH_with_traps as parent of different subtypes.

"""
abstract type model_SRH_with_traps  end

abstract type model_SRH_traps_transient <: model_SRH_with_traps    end


##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for model type which indicates, if we consider stationary
or transient problem.

"""
abstract type model_type end


"""
$(TYPEDEF)
Abstract type for transient simulations.

"""
abstract type model_transient <: model_type end


"""
$(TYPEDEF)
Abstract type for stationary simulations.

"""
abstract type model_stationary <: model_type end
##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for flux discretization model which is a parent of several subtypes.

"""
abstract type flux_approximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization.
Choose this one, when the Boltzmann statistics function is
chosen as statistics, check D. Scharfetter and H. Gummel, “Large-signal analysis of a silicon
Read diode oscillator”, IEEE Trans. Electr. Dev., vol. 16, pp. 64–77, 1969.

"""
abstract type scharfetter_gummel <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization for graded
effective density of states and/or graded band-edge energies. This means,
use this flux when at least one of these quantities
is assumed to be space-dependent.

"""
abstract type scharfetter_gummel_graded <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization, check  Z. Yu, and R. Dutton,
“SEDAN III – A one-dimensional device simulator”,
http://www-tcad.stanford.edu/tcad/programs/sedan3.html, 1988.

"""
abstract type excess_chemical_potential <: flux_approximation end

"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization
for graded effective density of states and/or graded band-edge
energies. This means, use this flux when at least one of these quantities
is assumed to be space-dependent.

"""
abstract type excess_chemical_potential_graded <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for diffusion enhanced flux discretization, check
M. Bessemoulin-Chatard, “A finite volume scheme for convection–diffusion equations with
nonlinear diffusion derived from the Scharfetter–Gummel scheme”, Numerische Mathematik,
vol. 121, pp. 637–670, 2012.

"""
abstract type diffusion_enhanced <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for generalized Scharfetter-Gummel flux discretization.
This flux approximation results in an implicit equation which needs to be
solved and is exact for all Blakemore type statistics functions with
abritary γ, check T. Koprucki and K. Gärtner. “Discretization scheme for drift-diffusion
equations with strong diffusion enhancement”. In: 12th International Conference on Numerical
Simulation of Optoelectronic Devices (NUSOD). 2012, pp. 103–104.

"""
abstract type generalized_sg <: flux_approximation end

##########################################################
##########################################################
"""
$(TYPEDEF)

Abstract type calculation_type which distinguishes between equilibrium and out
of equilibrium calculations.

"""
abstract type calculation_type end


"""
$(TYPEDEF)

Abstract type for equilibrium calculations.

"""
abstract type inEquilibrium <: calculation_type end


"""
$(TYPEDEF)

Abstract type for outOfEquilibrium calculations.

"""
abstract type outOfEquilibrium <: calculation_type end

##########################################################
##########################################################
"""
$(TYPEDEF)

Abstract type for generation model which is a parent of several subtypes.

"""
abstract type generation_model end


"""
$(TYPEDEF)

Abstract type for uniform generation.

"""
abstract type generation_uniform <: generation_model end


"""
$(TYPEDEF)

Abstract type for Beer-Lambert generation. Note that this type is implemented, but
not well tested yet.

"""
abstract type generation_beer_lambert <: generation_model end


"""
$(TYPEDEF)

Abstract type for no generation model.

"""
abstract type generation_none <: generation_model end

##########################################################
##########################################################

"""
$(TYPEDEF)
Abstract type for scan protocol type.

"""
abstract type scan_protocol_type end


"""
$(TYPEDEF)
Abstract type for linear scan protocol.

"""
abstract type linearScanProtocol <: scan_protocol_type end