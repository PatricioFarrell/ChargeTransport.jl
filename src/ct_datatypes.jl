
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
Abstract type for ohmic contacts as boundary model.

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
is part of boundary model. Subtypes are given below.

"""
abstract type interface_model <: boundary_model end


abstract type interface_model_tangential_flux                           <: interface_model end

abstract type interface_model_surface_recombination_and_tangential_flux <: interface_model end

"""
$(TYPEDEF)
Abstract type for no interface model.

"""
abstract type interface_model_none <: interface_model end


"""
$(TYPEDEF)
Abstract type for an interface model where discontinuous 
quasi Fermi potentials are needed.

"""
abstract type interface_model_discont_qF <: interface_model end


"""
$(TYPEDEF)
Abstract type for surface recombination mechanisms.

"""
abstract type interface_model_surface_recombination <: interface_model end


"""
$(TYPEDEF)
Abstract type for present ion charges at interfaces.

"""
abstract type interface_model_ion_charge <: interface_model end

##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for SRH bulk recombination model

    !!! compat  
    This one will be removed in future versions.

"""
abstract type abstract_model_SRH end

"""
$(TYPEDEF)
model_SRH as parent of several different subtypes.

"""
abstract type model_SRH                          <: abstract_model_SRH      end


"""
$(TYPEDEF)
model_SRH_without_traps as parent of several different subtypes.

"""
abstract type model_SRH_without_traps            <: model_SRH               end

abstract type model_SRH_stationary               <: model_SRH_without_traps end
abstract type model_SRH_off                      <: model_SRH_without_traps end

"""
$(TYPEDEF)
model_SRH_with_traps as parent of several different subtypes.

"""
abstract type model_SRH_with_traps               <: model_SRH               end

abstract type model_SRH_traps_transient          <: model_SRH_with_traps    end



"""
$(TYPEDEF)
This Datatype will be deleted soon.

"""
abstract type model_SRH_2species_present_trap_dens <: abstract_model_SRH end

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
Abstract type for flux discretization model.

"""
abstract type flux_approximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization.
Choose this one, when the Boltzmann statistics function is
chosen as statistics.

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
Abstract type for excess chemical potential flux discretization.

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
Abstract type for diffusion enhanced flux discretization.

"""
abstract type diffusion_enhanced <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for generalized Scharfetter-Gummel flux discretization.
This flux approximation results in an implicit equation which needs to be
solved and is exact for all Blakemore type statistics functions with
abritary γ.

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

Abstract type for out of equilibrium calculations.

"""
abstract type outOfEquilibrium <: calculation_type end

##########################################################
##########################################################
"""
$(TYPEDEF)

Abstract type for generation model.

"""
abstract type generation_model end


"""
$(TYPEDEF)

Abstract type for uniform generation.

"""
abstract type generation_uniform <: generation_model end


"""
$(TYPEDEF)

Abstract type for Beer-Lambert generation.

"""
abstract type generation_beer_lambert <: generation_model end


"""
$(TYPEDEF)

Abstract type for no generation model.

"""
abstract type generation_none <: generation_model end