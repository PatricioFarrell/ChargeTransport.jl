##########################################################
##########################################################
"""

$(TYPEDEF)
Abstract type for boundary model. Subtypes are OhmicContact and SchottkyContact.

"""
abstract type BoundaryModel   end

############    outer boundary conditions     ############
"""
Abstract type for ohmic contacts as outer boundary model.

"""
abstract type OhmicContact <: BoundaryModel  end


"""
Abstract type for schottky contacts as boundary model.

"""
abstract type SchottkyContact <: BoundaryModel end

############    inner boundary conditions     ############
"""
$(TYPEDEF)
Abstract type for interface model with several subtypes.

"""
abstract type InterfaceModel end

abstract type InterfaceModelTangentialFlux               <: InterfaceModel end

abstract type InterfaceModelSurfaceRecoAndTangentialFlux <: InterfaceModel end

"""
$(TYPEDEF)
Abstract type for no interface model.

"""
abstract type InterfaceModelNone <: InterfaceModel end

#
#$(TYPEDEF)
#Abstract type for an interface model where discontinuous
#quasi Fermi potentials are needed.
#
abstract type InterfaceModelDiscontqF <: InterfaceModel end


"""
$(TYPEDEF)
Abstract type for surface recombination mechanisms.

"""
abstract type InterfaceModelSurfaceReco <: InterfaceModel end


# """
# $(TYPEDEF)
# Abstract type for present ion charges at interfaces.

# """
abstract type InterfaceModelIonCharge <: InterfaceModel end

##########################################################
##########################################################

"""
$(TYPEDEF)
SRHWithoutTraps as parent of different subtypes.

"""
abstract type SRHWithoutTraps end

abstract type SRHStationary     <: SRHWithoutTraps end
abstract type SRHOff            <: SRHWithoutTraps end

"""
$(TYPEDEF)
SRHWithTraps as parent of different subtypes.

"""
abstract type SRHWithTraps  end

abstract type SRHTrapsTransient <: SRHWithTraps end


##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for model type which indicates, if we consider stationary
or transient problem.

"""
abstract type ModelType end


"""
$(TYPEDEF)
Abstract type for transient simulations.

"""
abstract type Transient <: ModelType end


"""
$(TYPEDEF)
Abstract type for stationary simulations.

"""
abstract type Stationary <: ModelType end
##########################################################
##########################################################
"""
$(TYPEDEF)
Abstract type for flux discretization model which is a parent of several subtypes.

"""
abstract type FluxApproximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization.
Choose this one, when the Boltzmann statistics function is
chosen as statistics, check D. Scharfetter and H. Gummel, “Large-signal analysis of a silicon
Read diode oscillator”, IEEE Trans. Electr. Dev., vol. 16, pp. 64–77, 1969.

"""
abstract type ScharfetterGummel <: FluxApproximation end


"""
$(TYPEDEF)
Abstract type for Scharfetter-Gummel flux discretization for graded
effective density of states and/or graded band-edge energies. This means,
use this flux when at least one of these quantities
is assumed to be space-dependent.

"""
abstract type ScharfetterGummelGraded <: FluxApproximation end


"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization, check  Z. Yu, and R. Dutton,
“SEDAN III – A one-dimensional device simulator”,
http://www-tcad.stanford.edu/tcad/programs/sedan3.html, 1988.

"""
abstract type ExcessChemicalPotential <: FluxApproximation end

"""
$(TYPEDEF)
Abstract type for excess chemical potential flux discretization
for graded effective density of states and/or graded band-edge
energies. This means, use this flux when at least one of these quantities
is assumed to be space-dependent.

"""
abstract type ExcessChemicalPotentialGraded <: FluxApproximation end


"""
$(TYPEDEF)
Abstract type for diffusion enhanced flux discretization, check
M. Bessemoulin-Chatard, “A finite volume scheme for convection–diffusion equations with
nonlinear diffusion derived from the Scharfetter–Gummel scheme”, Numerische Mathematik,
vol. 121, pp. 637–670, 2012.

"""
abstract type DiffusionEnhanced <: FluxApproximation end


"""
$(TYPEDEF)
Abstract type for generalized Scharfetter-Gummel flux discretization.
This flux approximation results in an implicit equation which needs to be
solved and is exact for all Blakemore type statistics functions with
abritary γ, check T. Koprucki and K. Gärtner. “Discretization scheme for drift-diffusion
equations with strong diffusion enhancement”. In: 12th International Conference on Numerical
Simulation of Optoelectronic Devices (NUSOD). 2012, pp. 103–104.

"""
abstract type GeneralizedSG <: FluxApproximation end

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
##########################################################
"""
$(TYPEDEF)

Abstract type for generation model which is a parent of several subtypes.

"""
abstract type GenerationModel end


"""
$(TYPEDEF)

Abstract type for uniform generation.

"""
abstract type GenerationUniform <: GenerationModel end


"""
$(TYPEDEF)

Abstract type for Beer-Lambert generation. Note that this type is implemented, but
not well tested yet.

"""
abstract type GenerationBeerLambert <: GenerationModel end


"""
$(TYPEDEF)

Abstract type for no generation model.

"""
abstract type GenerationNone <: GenerationModel end

##########################################################
##########################################################

"""
$(TYPEDEF)
Abstract type for scan protocol type.

"""
abstract type ScanProtocolType end


"""
$(TYPEDEF)
Abstract type for linear scan protocol.

"""
abstract type LinearScanProtocol <: ScanProtocolType end