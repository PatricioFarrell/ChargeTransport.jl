"""
Type of statistics functions.

"""
const StandardFuncSet = Union{typeof(Boltzmann), typeof(Blakemore), typeof(FermiDiracMinusOne),
                              typeof(FermiDiracOneHalfBednarczyk), typeof(FermiDiracOneHalfTeSCA)}

##########################################################
"""
Type of charge carriers and the electric potential.

"""
const QType = Union{VoronoiFVM.ContinuousQuantity{Int64}, VoronoiFVM.DiscontinuousQuantity{Int64}, Int64}

##########################################################
##########################################################

############    outer boundary conditions     ############
"""
Abstract type for ohmic contacts as outer boundary model.

"""
abstract type OhmicContact  end


"""
Abstract type for schottky contacts as boundary model.

"""
abstract type SchottkyContact end

############    inner boundary conditions     ############
"""
$(TYPEDEF)
Abstract type for no interface model.

"""
abstract type InterfaceModelNone end

#
#$(TYPEDEF)
#Abstract type for an interface model where discontinuous
#quasi Fermi potentials are needed.
#
abstract type InterfaceModelDiscontqF end

"""
$(TYPEDEF)
Abstract type for surface recombination mechanisms.

"""
abstract type InterfaceModelSurfaceReco end

abstract type InterfaceModelTangentialFlux end

abstract type InterfaceModelSurfaceRecoAndTangentialFlux end

# """
# $(TYPEDEF)
# Abstract type for present ion charges at interfaces.

# """
abstract type InterfaceModelIonCharge end

##########################################################
"""
Possible types of outer boundary model.

"""
const OuterBoundaryModelType = Union{Type{OhmicContact}, Type{SchottkyContact}}


"""
Possible Types of interface model (interior boundary conditions).

"""
const InterfaceModelType = Union{Type{InterfaceModelNone}, Type{InterfaceModelSurfaceReco},
                             Type{InterfaceModelDiscontqF}, Type{InterfaceModelTangentialFlux},
                             Type{InterfaceModelSurfaceRecoAndTangentialFlux}, Type{InterfaceModelIonCharge}}

"""
Possible types of boundary models.

"""
const BoundaryModelType  = Union{OuterBoundaryModelType, InterfaceModelType}

##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for transient simulations.

"""
abstract type Transient end


"""
$(TYPEDEF)
Abstract type for stationary simulations.

"""
abstract type Stationary end

##########################################################
"""
Possible types which indicate, if we consider stationary or transient problem.

"""
const ModelType = Union{Type{Transient}, Type{Stationary}}

##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization. Choose this one, when the Boltzmann
statistics function is chosen as statistics, check
D. Scharfetter and H. Gummel, “Large-signal analysis of a silicon Read diode oscillator”, IEEE Trans. Electr. Dev., vol. 16, pp. 64–77, 1969.

"""
abstract type ScharfetterGummel end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization for graded effective density of
states and/or graded band-edge energies. This means, use this flux when at least one of these
parameters is assumed to be space-dependent.

"""
abstract type ScharfetterGummelGraded end

"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization, check  Z. Yu, and R. Dutton,
“SEDAN III – A one-dimensional device simulator”,
http://www-tcad.stanford.edu/tcad/programs/sedan3.html, 1988.

"""
abstract type ExcessChemicalPotential end

"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization for graded effective density
of states and/or graded band-edge energies. This means, use this flux when at least one of
these parameters is assumed to be space-dependent.

"""
abstract type ExcessChemicalPotentialGraded end

"""
$(TYPEDEF)
Abstract type for diffusion enhanced flux discretization, check
M. Bessemoulin-Chatard, “A finite volume scheme for convection–diffusion equations with
nonlinear diffusion derived from the Scharfetter–Gummel scheme”, Numerische Mathematik,
vol. 121, pp. 637–670, 2012.

"""
abstract type DiffusionEnhanced end

"""
$(TYPEDEF)
Abstract type for generalized Scharfetter-Gummel flux discretization. This flux approximation
results in an implicit equation which needs to be solved and is exact for all Blakemore type
statistics functions with abritary γ, check T. Koprucki and K. Gärtner.
“Discretization scheme for drift-diffusion equations with strong diffusion enhancement”.
In: 12th International Conference on Numerical Simulation of Optoelectronic Devices (NUSOD). 2012, pp. 103–104.

"""
abstract type GeneralizedSG end

##########################################################
"""
Possible types of flux discretization schemes.

"""
const FluxApproximationType = Union{Type{ScharfetterGummel}, Type{ExcessChemicalPotential},
                                 Type{DiffusionEnhanced}, Type{GeneralizedSG},
                                 Type{ScharfetterGummelGraded}, Type{ExcessChemicalPotentialGraded}}

##########################################################
##########################################################
"""
$(TYPEDEF)

Abstract type for equilibrium calculations.

"""
abstract type InEquilibrium end


"""
$(TYPEDEF)

Abstract type for out of equilibrium calculations.

"""
abstract type OutOfEquilibrium end

##########################################################
"""
Possible types for calculation type.

"""
const CalculationType = Union{Type{InEquilibrium}, Type{OutOfEquilibrium}}

##########################################################
##########################################################
abstract type AbstractModelSRH end

abstract type ModelSRH <: AbstractModelSRH      end
abstract type SRH2SpeciesPresentTrapDens <: AbstractModelSRH end
##############################################

abstract type SRHStationary <:ModelSRH end
abstract type SRHOff       <:ModelSRH  end

abstract type SRHTrapsTransient end

##########################################################
"""
Possible type for SRH recombination without traps.

"""
const SRHWithoutTrapsType = Union{Type{SRHStationary}, Type{SRHOff}}

"""
Possible types for SRH recombination without traps.

"""
const SRHWithTrapsType = Type{SRHTrapsTransient}

const SRHModelType = Union{SRHWithoutTrapsType, SRHWithTrapsType}

##########################################################
##########################################################

"""
$(TYPEDEF)

Abstract type for uniform generation.

"""
abstract type GenerationUniform end


"""
$(TYPEDEF)

Abstract type for Beer-Lambert generation. Note that this type is implemented, but
not well tested yet.

"""
abstract type GenerationBeerLambert end


"""
$(TYPEDEF)

Abstract type for no generation model.

"""
abstract type GenerationNone end

##########################################################
"""
Possible types for generation model.

"""
const GenerationModelType = Union{Type{GenerationUniform}, Type{GenerationBeerLambert}, Type{GenerationNone}}


##########################################################