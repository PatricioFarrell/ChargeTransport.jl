##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for boundary model.

"""
abstract type boundary_model   end 

############    outer boundary conditions     ############
"""
$(TYPEDEF)
Abstract type for ohmic contacts as boundary model.

"""
abstract type ohmic_contact <: boundary_model  end


"""
$(TYPEDEF)
Abstract type for schottky contacts as boundary model.

"""
abstract type schottky_contact <: boundary_model end

############    inner boundary conditions     ############
"""
$(TYPEDEF)
Abstract type for interface model which
is part of boundary model.

"""
abstract type interface_model <: boundary_model end


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


"""
$(TYPEDEF)
Abstract type for distinction of ion charge interface model
between left and right of active layer. Distinguishing is necessary
due to sign of electrochemical reaction.

"""
abstract type interface_model_ion_charge_left <: interface_model_ion_charge end


"""
$(TYPEDEF)
Abstract type for distinction of ion charge interface model
between left and right of active layer. Distinguishing is necessary
due to sign of electrochemical reaction.

"""
abstract type interface_model_ion_charge_right <: interface_model_ion_charge end

##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for bulk recombination model.

"""
abstract type bulk_recombination_model end


"""
$(TYPEDEF)
Abstract type for no bulk recombination model.

"""
abstract type bulk_recomb_model_none <: bulk_recombination_model end 


"""
$(TYPEDEF)
Abstract type for trap assisted bulk recombination model, i.e.
only Schockley-Read-Hall recombination is used.

"""
abstract type bulk_recomb_model_trap_assisted <: bulk_recombination_model end


"""
$(TYPEDEF)
Abstract type for only radiative recombination model.

"""
abstract type bulk_recomb_model_radiative <: bulk_recombination_model end


"""
$(TYPEDEF)
Abstract type for full bulk recombination model.
Currently, Schockley-Read-Hall, radiative and Auger are implemented.

"""
abstract type bulk_recomb_model_full <: bulk_recombination_model end

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
abstract type ScharfetterGummel <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization for graded
effective density of states and/or graded band-edge energies, i.e.
when this two quantities are assumed to be space-dependent.

"""
abstract type ScharfetterGummel_Graded <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization.

"""
abstract type excessChemicalPotential <: flux_approximation end

"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization
for graded effective density of states and/or graded band-edge 
energies, i.e. when this two quantities are assumed to be space-dependent.

"""
abstract type excessChemicalPotential_Graded <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for diffusion enhanced flux discretization.

"""
abstract type diffusionEnhanced <: flux_approximation end


"""
$(TYPEDEF)
Abstract type for generalized Scharfetter-Gummel flux discretization.
This flux approximation results in an implicit equation which needs to be
solved and is exact for all Blakemore type statistics functions with
abritary γ.

"""
abstract type generalized_SG <: flux_approximation end

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