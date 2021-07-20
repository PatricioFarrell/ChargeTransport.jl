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
abstract type bulk_recombination_none <: bulk_recombination_model end 


"""
$(TYPEDEF)
Abstract type for trap assisted bulk recombination model, i.e.
only Schockley-Read-Hall recombination is used.

"""
abstract type bulk_recombination_trap_assisted <: bulk_recombination_model end


"""
$(TYPEDEF)
Abstract type for only radiative recombination model.

"""
abstract type bulk_recombination_radiative <: bulk_recombination_model end


"""
$(TYPEDEF)
Abstract type for full bulk recombination model.
Currently, Schockley-Read-Hall, radiative and Auger are implemented.

"""
abstract type bulk_recombination_full <: bulk_recombination_model end

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
##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for interior nodes.
"""
function etaFunction(u, node::VoronoiFVM.Node, data, icc::Int64, ipsi::Int64)
    params      = data.params
    paramsnodal = data.paramsnodal

    E  = params.bandEdgeEnergy[icc, node.region] + paramsnodal.bandEdgeEnergy[icc, node.index]
    
    return params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for boundary nodes.
"""
function etaFunction(u, bnode::VoronoiFVM.BNode, data, icc::Int64, ipsi::Int64) # bnode.index refers to index in overall mesh
    params      = data.params
    paramsnodal = data.paramsnodal

    E  = params.bBandEdgeEnergy[icc, bnode.region] + paramsnodal.bandEdgeEnergy[icc, bnode.index]
    
    return params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for edges.
"""

function etaFunction(u, edge::VoronoiFVM.Edge, data, icc::Int64, ipsi::Int64, nodeEdge)

    params      = data.params
    paramsnodal = data.paramsnodal
    
    E  = params.bandEdgeEnergy[icc, edge.region] + paramsnodal.bandEdgeEnergy[icc, nodeEdge]

    return params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for given ``\\varphi_\\alpha``
and ``\\psi``

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q ).``

The parameters ``E_\\alpha`` and ``z_\\alpha`` are given as vectors.
This function may be used to compute the charge density, i.e. the
right-hand side of the Poisson equation.   
"""
function etaFunction(psi, phi, UT, E::Array, z::Array)
    z ./ UT .* ( (phi - psi) .+ E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for floats.
"""
function etaFunction(u, data, node, region, icc::Int64, ipsi::Int64, in_region::Bool)

    params      = data.params
    paramsnodal = data.paramsnodal

    if in_region == true
        E  = params.bandEdgeEnergy[icc, region] + paramsnodal.bandEdgeEnergy[icc, node]
    elseif in_region == false
        E  = params.bBandEdgeEnergy[icc, region] + paramsnodal.bandEdgeEnergy[icc, node]
    end

    return params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
end

##########################################################
##########################################################

function emptyFunction()
end


"""
$(TYPEDSIGNATURES)
Master breaction! function. This is the function which enters VoronoiFVM and hands over
for each boundary the chosen boundary model.

"""
breaction!(f, u, bnode, data) =  breaction!(f, u, bnode, data, data.boundary_type[bnode.region])



"""
$(TYPEDSIGNATURES)

This breaction! function is chosen when no interface model is chosen.

"""
breaction!(f, u, bnode, data, ::Type{interface_model_none}) = emptyFunction()



"""
$(TYPEDSIGNATURES)

Creates ohmic boundary conditions via a penalty approach with penalty parameter ``\\delta``.
For example, the right-hand side for the electrostatic potential ``\\psi`` is implemented as

``f[\\psi]  = -q/\\delta   ( (p - N_a) - (n - N_d) )``,

assuming a bipolar semiconductor. In general, we have for some given charge number ``z_\\alpha``

``f[\\psi] =  -q/\\delta  \\sum_\\alpha{ z_\\alpha  (n_\\alpha - C_\\alpha) },``

where ``C_\\alpha`` corresponds to some doping w.r.t. the species ``\\alpha``.

The boundary conditions for the charge carriers are set in the main file. Hence,

``f[n_\\alpha] = 0```

for all charge carriers ``\\alpha``.
"""
function breaction!(f, u, bnode, data, ::Type{ohmic_contact})

    params      = data.params
    paramsnodal = data.paramsnodal 

     # parameters
     α          = 1.0e-10                      # tiny penalty value
     ipsi       = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1  # final index for electrostatic potential
 
 
    for icc = 1:params.numberOfCarriers
 
        eta     = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )                    # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.bDensityOfStates[icc, bnode.region] + paramsnodal.densityOfStates[icc, bnode.index]) * data.F[icc](eta)  # add charge carrier
 
        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0
 
    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]

    f[ipsi] = - data.λ1 * 1 / α *  q * f[ipsi]

    return f

end


"""
$(TYPEDSIGNATURES)

This breaction! function is chosen  when we assume 
ion charges to be present at inner interfaces for the left inner boundary.

"""
function breaction!(f, u, bnode, data, ::Type{interface_model_ion_charge_left})  

    params            = data.params
    paramsnodal       = data.paramsnodal 

    iphia             = 3
    iphiaj1, iphiaj2  = 4:5
    ipsi              = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1

    E1                = params.bBandEdgeEnergy[iphia, bnode.region]  + paramsnodal.bandEdgeEnergy[iphia, bnode.index]
    DOS1              = params.bDensityOfStates[iphia, bnode.region] + paramsnodal.densityOfStates[iphia, bnode.index]
    C01               = params.bDoping[iphia, bnode.region]  

    β                 = 0.5     # can be between 0 and 1 
    κ                 = 1       # either 0 or 1
    r0                = params.r0

    etaInterfaceAnion = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaj1] - u[ipsi]) + E1 / q )
    
    f[ipsi]           =  - data.λ1 * q  * ( params.chargeNumbers[iphia] * DOS1^(2/3) * data.F[iphia](etaInterfaceAnion) - C01^(2/3) ) # (1.4.5) @ left inner boundary 
    
    # DA: das kann noch besser gemacht werden ....
    if data.calculation_type == inEquilibrium
        f[iphia]       = u[iphia]
        f[iphiaj1]     = u[iphiaj1]
     else
    
        f[iphia]       =   data.λ3 * q * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj1, ipsi, β, κ, DOS1, E1) ) # (1.4.8) @ left inner boundary 
        f[iphiaj1]     = - data.λ3 * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj1, ipsi, β, κ, DOS1, E1) ) # (1.4.7) @ left inner boundary (right-hand side of equation)
    end

    return f
end

"""
$(TYPEDSIGNATURES)

This breaction! function is chosen  when we assume 
ion charges to be present at inner interfaces for the right inner boundary.

"""
function breaction!(f, u, bnode, data, ::Type{interface_model_ion_charge_right})

    params            = data.params
    paramsnodal       = data.paramsnodal 

    iphia             = 3
    iphiaj1, iphiaj2  = 4:5
    ipsi              = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1

    E2                = params.bBandEdgeEnergy[iphia, bnode.region]  + paramsnodal.bandEdgeEnergy[iphia, bnode.index]
    DOS2              = params.bDensityOfStates[iphia, bnode.region] + paramsnodal.densityOfStates[iphia, bnode.index]
    C02               = params.bDoping[iphia, bnode.region]         

    β                 = 0.5     # can be between 0 and 1 
    κ                 = 1       # either 0 or 1
    r0                = params.r0

    etaInterfaceAnion = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaj2] - u[ipsi]) + E2 / q )
    
    f[ipsi]           =  - data.λ1 * q * ( params.chargeNumbers[iphia] * DOS2^(2/3) * data.F[iphia](etaInterfaceAnion) - C02^(2/3) ) # (1.4.5) @ rigth inner boundary 


    if data.calculation_type == inEquilibrium
        f[iphia]       = u[iphia]
        f[iphiaj2]     = u[iphiaj2]
    else

        f[iphia]       = - data.λ3 *  q * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj2, ipsi, β, κ, DOS2, E2) ) # (1.4.8) @ right inner boundary 
        f[iphiaj2]     = - data.λ3 * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj2, ipsi, β, κ, DOS2, E2) ) # (1.4.7) @ right inner boundary (right-hand side of equation)
    end

    return f
    
end


"""
$(TYPEDSIGNATURES)

Electrochemical reaction between interface and bulk ionic species.
This function enters in the internal boundary reaction
in case of an ion charge interface model.

"""
function electrochemicalReaction(data, u, iphia, ipsi, iphiaJunction, ipsiJunction, β, κ, DOS, E) # (1.4.9)

    params             = data.params
 

    etaExp             = params.chargeNumbers[iphia] / params.UT * ( (u[iphia] - u[iphiaJunction]) + E / q ) 
    expTerm            =  exp( β * etaExp ) - exp( (β - 1) * etaExp)

    etaInterfaceAnion  = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaJunction] - u[ipsiJunction]) + E / q )
    etaAnion           = params.chargeNumbers[iphia] / params.UT * ( (u[iphia] - u[ipsi]) + E / q )

    densFactor         = ( (DOS^(2/3) * data.F[iphia](etaInterfaceAnion) )^(1/2) * (DOS * data.F[iphia](etaAnion) )^(- 1/2) )^κ

    return densFactor * expTerm

end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bstorage! function. This is the function which enters VoronoiFVM and hands over
for each boundary the time-dependent part of the chosen boundary model.

"""
bstorage!(f, u, bnode, data) = bstorage!(f, u, bnode, data, data.boundary_type[bnode.region])


"""
$(TYPEDSIGNATURES)
No bstorage! is used, if no interface model is chosen.

"""
bstorage!(f, u, bnode, data, ::Type{interface_model_none}) = emptyFunction()

"""
$(TYPEDSIGNATURES)
No bstorage! is used, if surface recombination model is chosen.

"""
bstorage!(f, u, bnode, data, ::Type{interface_model_surface_recombination}) = emptyFunction()

"""
$(TYPEDSIGNATURES)
No bstorage! is used, if an ohmic contact model is chosen.

"""
bstorage!(f, u, bnode, data, ::Type{ohmic_contact}) = emptyFunction()



"""
$(TYPEDSIGNATURES)
Time-dependent part in case of present ion charges at inner interfaces (left).

"""
function bstorage!(f, u, bnode, data, ::Type{interface_model_ion_charge_left})
# DA: I guess, here is no sign included so we do not even need to distinguish between left and right.
# But need to fix use of interface charges (done by indices)!

    params            = data.params
    paramsnodal        = data.paramsnodal

    iphia             = 3
    iphiaj1, iphiaj2  = 4:5
    ipsi              = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1        # final index for electrostatic potential

    E1                = params.bBandEdgeEnergy[iphia, bnode.region]  + paramsnodal.bandEdgeEnergy[iphia, bnode.index]
    DOS1              = params.bDensityOfStates[iphia, bnode.region] + paramsnodal.densityOfStates[iphia, bnode.index]

    # (1.4.7) @ left inner boundary (left-hand side of equation)
    etaInterfaceAnion = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaj1] - u[ipsi]) + E1 / q ) 
    f[iphiaj1]        = params.chargeNumbers[iphia] * DOS1^(2/3) * data.F[iphia](etaInterfaceAnion) 

    f[ipsi]           = 0.0

    return f

end


"""
$(TYPEDSIGNATURES)
Time-dependent part in case of present ion charges at inner interfaces (right).

"""
function bstorage!(f, u, bnode, data, ::Type{interface_model_ion_charge_right})

    params            = data.params
    paramsnodal        = data.paramsnodal

    iphia             = 3
    iphiaj1, iphiaj2  = 4:5
    ipsi              = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1        # final index for electrostatic potential

    
    E2                = params.bBandEdgeEnergy[iphia, bnode.region]  + paramsnodal.bandEdgeEnergy[iphia, bnode.index]
    DOS2              = params.bDensityOfStates[iphia, bnode.region] + paramsnodal.densityOfStates[iphia, bnode.index]

    # (1.4.7) @ right inner boundary (left-hand side of equation)
    etaInterfaceAnion = params.chargeNumbers[iphia] / params.UT * ( (u[iphiaj2] - u[ipsi]) + E2 / q )
    f[iphiaj2]        = params.chargeNumbers[iphia] *  DOS2^(2/3) * data.F[iphia](etaInterfaceAnion)

    f[ipsi]           = 0.0

    return f

end


##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master reaction! function. This is the function which enters VoronoiFVM and hands over
reaction terms for concrete calculation type and bulk recombination model.

"""
reaction!(f, u, node, data) = reaction!(f, u, node, data, data.calculation_type)

"""
$(TYPEDSIGNATURES)
Reaction in case of equilibrium, i.e. no generation and recombination is considered.

"""
function reaction!(f, u, node, data, ::Type{inEquilibrium})

    params      = data.params
    paramsnodal = data.paramsnodal

    # indices
    ipsi                  = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1            # final index for electrostatic potential
    ireg                  = node.region
    inode                 = node.index

    # rhs of NLP (charge density)
    for icc = 1:params.numberOfCarriers
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, ireg] )   # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, inode]) * data.F[icc](eta)   # add charge carrier

    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[inode]

    f[ipsi]                     = - data.λ1 * q * f[ipsi]

    for icc = 1:params.numberOfCarriers
        f[icc] = u[icc] - 0.0
    end

    return f

end

recombination_kernel(data, ireg, iphin::Int64, iphip::Int64, n, p, ::Type{bulk_recombination_none}) = 0.0


function recombination_kernel(data, ireg, iphin::Int64, iphip::Int64,  n, p, ::Type{bulk_recombination_radiative})

    params = data.params

    return params.recombinationRadiative[ireg]

end


function recombination_kernel(data, ireg, iphin::Int64, iphip::Int64, n, p,::Type{bulk_recombination_trap_assisted})

    params = data.params

    kernelSRH = 1.0 / (  params.recombinationSRHLifetime[iphip, ireg] * (n + params.recombinationSRHTrapDensity[iphin, ireg]) + params.recombinationSRHLifetime[iphin, ireg] * (p + params.recombinationSRHTrapDensity[iphip, ireg]) )

    return kernelSRH

end


function recombination_kernel(data, ireg, iphin::Int64, iphip::Int64, n, p, ::Type{bulk_recombination_full})

    params = data.params

    # radiative recombination
    kernelRadiative = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recombination_radiative)

    # SRH recombination
    kernelSRH       = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recombination_trap_assisted)

    # Auger recombination
    kernelAuger     = (params.recombinationAuger[iphin, ireg] * n + params.recombinationAuger[iphip, ireg] * p)


    return kernelRadiative + kernelAuger + kernelSRH

end

"""
$(TYPEDSIGNATURES)

Sets up the right-hand sides. Assuming a bipolar semiconductor
the right-hand side for the electrostatic potential becomes

  ``f[ψ]  = - q ((p - N_a) - (n - N_d) ) = - q  \\sum  n_\\alpha  (n_\\alpha - C_\\alpha) ``

for some doping ``C_\\alpha`` w.r.t. to the species ``\\alpha``.
The right-hand sides for the charge carriers read as

``f[n_\\alpha] =  - z_\\alpha  q (G -  R) ``

for all charge carriers ``n_\\alpha``.
The recombination includes radiative, Auger and Shockley-Read-Hall
recombination. For latter recombination process the stationary simplification is implemented.

The recombination is only implemented for electron and holes and assumes
that the electron index is 1 and the hole index is 2. 

"""
function reaction!(f, u, node, data, ::Type{outOfEquilibrium})

    params      = data.params
    paramsnodal = data.paramsnodal

    # indices
    iphin                 = 1
    iphip                 = 2
    ipsi                  = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1            # final index for electrostatic potential
    ireg                  = node.region
    inode                 = node.index

    # rhs of NLP (charge density)
    for icc = 1:params.numberOfCarriers
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[inode]

    # rhs of continuity equations for electron and holes (bipolar reaction)
    n                     = compute_densities!(u, data, inode, ireg, iphin, ipsi, true)  # true for interior region
    p                     = compute_densities!(u, data, inode, ireg, iphip, ipsi, true) 
    exponentialTerm       = exp((q *u[iphin] - q  * u[iphip]) / (kB*params.temperature))
    excessCarrierDensTerm = n*p * (1.0 - exponentialTerm)

    for icc in [iphin, iphip] 

        # gives you the recombination kernel based on choice of user
        kernel = recombination_kernel(data, ireg, iphin, iphip, n, p, data.bulk_recombination_model)
                
        f[icc]          = q * params.chargeNumbers[icc] *  kernel *  excessCarrierDensTerm  - q * params.chargeNumbers[icc] * generation(data, ireg,  node.coord[node.index], data.generation_model)
    end

    # vorsicht bei diesem part! hat auswirkungen auf tier 0
    for icc in iphip+1:params.numberOfCarriers
        f[icc]              = u[icc] - 0.0
    end

    
    f[ipsi]                     = - q * data.λ1 * f[ipsi]

    return f

end


"""
$(SIGNATURES)

Compute trap densities for a given trap energy.
[Currently, only done for the Boltzmann statistics and for region dependent parameters.]

"""
function trap_density!(icc, ireg, data, Et) 
    params      = data.params

    params.densityOfStates[icc, ireg] * exp( params.chargeNumbers[icc] * (params.bandEdgeEnergy[icc, ireg] - Et) / (kB * params.temperature)) # need to subtract Eref
end


"""
$(TYPEDSIGNATURES)

The generation rate ``G``, which occurs in the right-hand side of the
continuity equations with a uniform generation rate.
"""
function generation(data, ireg, node, ::Type{generation_uniform}) # only works in 1D till now; adjust node, when multidimensions

    params = data.params

    return data.λ2 * params.generationUniform[ireg]
end


"""
$(TYPEDSIGNATURES)

The generation rate ``G``, which occurs in the right-hand side of the
continuity equations obeying the Beer-Lambert law.
"""
function generation(data, ireg, node, ::Type{generation_beer_lambert}) # only works in 1D till now; adjust node, when multidimensions

    params = data.params

    return data.λ2 * params.generationIncidentPhotonFlux[ireg] * params.generationAbsorption[ireg] * exp( - params.generationAbsorption[ireg] * node )

end

generation(data, ireg, node, ::Type{generation_none}) = 0.0

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master storage! function. This is the function which enters VoronoiFVM and hands over
a storage term, if we consider transient problem.

"""
storage!(f, u, node, data) = storage!(f, u, node, data, data.model_type)

storage!(f, u, node, data, ::Type{model_stationary})  = emptyFunction()


"""
$(TYPEDSIGNATURES)

The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have 

``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha`` 

and for the electrostatic potential
``f[ψ] = 0``.

"""
function storage!(f, u, node, data, ::Type{model_transient})

    params      = data.params
    paramsnodal = data.paramsnodal 

    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    
    for icc = 1:params.numberOfCarriers

        eta    = etaFunction(u, node, data, icc, ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = q * params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)
    end

    f[ipsi] =  0.0

    return f
end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Flux discretization scheme is chosen in two steps. First, we need
to see, if we are in or out of equilibrium. If, inEquilibrium, then
no flux is passed. If out of equilibrium, we choose the flux approximation
which the user chose.

"""

flux!(f, u, edge, data) = flux!(f, u, edge, data, data.calculation_type)

function flux!(f, u, edge, data, ::Type{inEquilibrium})

    params      = data.params
    paramsnodal  = data.paramsnodal

    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg        = edge.region

    dpsi        =   ul[ipsi] - uk[ipsi]
    f[ipsi]     = - (params.dielectricConstant[ireg] + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * ε0 * dpsi
    
end

flux!(f, u, edge, data, ::Type{outOfEquilibrium}) = flux!(f, u, edge, data, data.flux_approximation)


flux!(f, u, edge, data, ::Type{flux_approximation}) = emptyFunction()

"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{ScharfetterGummel})

    params      = data.params
    paramsnodal = data.paramsnodal

    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    dpsi        =   ul[ipsi] - uk[ipsi]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    for icc = 1:params.numberOfCarriers

        j0                 =  params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge)
        etal               = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm             = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT)
        f[icc]             = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
    
    end

    return f

end

"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme for 
possible space-dependent DOS and band-edge energies. For these parameters
the discretization scheme is modified. [insert continuous flux etc ...]

"""
function flux!(f, u, edge, data, ::Type{ScharfetterGummel_Graded})

    params      = data.params
    paramsnodal = data.paramsnodal

    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    dpsi        = ul[ipsi] - uk[ipsi]
    dpsiEps     = (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps
    
    for icc = 1:params.numberOfCarriers

        j0                 =  params.chargeNumbers[icc] * q * params.UT

        etak               = params.chargeNumbers[icc] / params.UT * ( (uk[icc] - uk[ipsi]) + (paramsnodal.bandEdgeEnergy[icc, nodek] + params.bandEdgeEnergy[icc, ireg]) / q )
        etal               = params.chargeNumbers[icc] / params.UT * ( (ul[icc] - ul[ipsi]) + (paramsnodal.bandEdgeEnergy[icc, nodel] + params.bandEdgeEnergy[icc, ireg]) / q )


        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobility           = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2
        densityOfStatesl   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, nodel])
        densityOfStatesk   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, nodek])
        
        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT ) 
        else
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc, nodek])) ) 
        end

        f[icc]             =  -j0  * mobility * ( bm  * densityOfStatesl * data.F[icc](etal) - bp *  densityOfStatesk * data.F[icc](etak) )

    end

    return f

end


"""
$(TYPEDSIGNATURES)

The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{excessChemicalPotential})

    params      = data.params
    paramsnodal = data.paramsnodal
    
    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    dpsi        =   ul[ipsi] - uk[ipsi]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    

    for icc = 1:params.numberOfCarriers

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u, edge, data, icc, ipsi, nodeEdge)
        etal               = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u, edge, data, icc, ipsi, nodeEdge)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q                  = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
    end

    return f

end


"""
$(TYPEDSIGNATURES)

The excess chemical potential flux scheme for 
possible space-dependent DOS and band-edge energies. For these parameters
the discretization scheme is modified. [insert continuous flux etc ...]

"""
function flux!(f, u, edge, data, ::Type{excessChemicalPotential_Graded})

    params      = data.params
    paramsnodal = data.paramsnodal

    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    dpsi        = ul[ipsi] - uk[ipsi]
    dpsiEps     = (params.dielectricConstant[ireg]  + paramsnodal.dielectricConstant[nodel]) * ul[ipsi] - (params.dielectricConstant[ireg] + paramsnodal.dielectricConstant[nodek]) * uk[ipsi]
    f[ipsi]     = - ε0 * dpsiEps
    
    
    for icc = 1:params.numberOfCarriers

        j0                 = params.chargeNumbers[icc] * q * params.UT
        etak               = params.chargeNumbers[icc] / params.UT * ( (uk[icc] - uk[ipsi]) + (paramsnodal.bandEdgeEnergy[icc, nodek] + params.bandEdgeEnergy[icc, ireg]) / q )
        etal               = params.chargeNumbers[icc] / params.UT * ( (ul[icc] - ul[ipsi]) + (paramsnodal.bandEdgeEnergy[icc, nodel] + params.bandEdgeEnergy[icc, ireg]) / q )

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobilityl          = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodel])
        mobilityk          = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodek])
        densityOfStatesl   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc,nodel])
        densityOfStatesk   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc,nodek])

        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            Q              = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )

        else
            Q              = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc,nodek])) )

        end

        bp, bm             = fbernoulli_pm(Q)
        f[icc]             = - j0  * ( bm  * mobilityl * densityOfStatesl * data.F[icc](etal) - bp * mobilityk * densityOfStatesk * data.F[icc](etak) )
    
    end

    return f

end


"""
$(TYPEDSIGNATURES)

The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is 
used for the regularization of the removable singularity. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{diffusionEnhanced})

    params      = data.params
    paramsnodal = data.paramsnodal

    tolReg      = 1.0e-13
    
    uk          = viewK(edge, u)
    ul          = viewL(edge, u)
    
    ipsi        = params.numberOfCarriers + params.numberOfInterfaceCarriers  + 1
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]
    
    dpsi        =   ul[ipsi] - uk[ipsi]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    for icc = 1:params.numberOfCarriers

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)
        etal               = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)

        if abs( (etal - etak)/(etak + etal) ) > tolReg
            g  = (etal - etak ) / ( log(data.F[icc](etal)) - log(data.F[icc](etak)) )
        else
            # regularization idea coming from Pietra-Jüngel scheme
            gk = exp(etak) / data.F[icc](etak)
            gl = exp(etal) / data.F[icc](etal)
            g  = 0.5 * ( gk + gl )
        end

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm             = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / (params.UT * g))
        f[icc]             = - j0 * g * (  bm * data.F[icc](etal) - bp * data.F[icc](etak))
    end

    return f

end

"""
$(TYPEDSIGNATURES)

The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation which arise
when considering the generalized Scharfetter-Gummel scheme in case of Blakemore statistics.
Hence, it should be exclusively worked with, when considering the Blakemore distribution.
This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{generalized_SG})

    params        = data.params
    paramsnodal   = data.paramsnodal

    max_iteration = 200          # for Newton solver
    it            = 0            # number of iterations (newton)
    damp          = 0.1          # damping factor
    
    uk            = viewK(edge, u)
    ul            = viewL(edge, u)
    nodel         = edge.node[2]
    nodek         = edge.node[1]
    
    ipsi          = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1
    ireg          = edge.region
    
    dpsi          =    ul[ipsi] - uk[ipsi]
    f[ipsi]       =  - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    for icc = 1:params.numberOfCarriers 

        j0                  = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        etak                = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)
        etal                = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)

        bandEdgeDifference  = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        # use Sedan flux as starting guess
        Q                   = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm              = fbernoulli_pm(Q)
        jInitial            = ( bm * data.F[icc](etal)  - bp * data.F[icc](etak))

        implicitEq(j::Real) = (fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q)) /params.UT + params.γ * j )[2] * exp(etal) - fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q) /params.UT) - params.γ * j)[1] * exp(etak)) - j

        delta               = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
        oldup               = 1.0
        while (it < max_iteration)
            Fval     = implicitEq(jInitial)
            dFval    = ForwardDiff.derivative(implicitEq, jInitial)

            if isnan(value(dFval)) || value(abs(dFval)) < delta
                @show value(jInitial), value(Fval), value(dFval)
                error("singular derivative in exact SG scheme")
            end
           
            update   = Fval / dFval
            jInitial = jInitial - damp * update

            if abs(update) < delta
                break
            end
            #@show abs(value(update)/oldup)
            oldup = value(update)

            it       = it + 1
            damp     = min(damp * 1.2, 1.0)
        end
        #@show  it
        f[icc]       = - j0 * jInitial
    end

    return f

end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
[Note that this way of implementation is not well tested yet. 
Will be coded different in future version !!!!]

Creates Schottky boundary conditions in a first attempt. For the electrostatic potential we assume 

``\\psi = \\psi_S + U, ``

where  ``\\psi_S`` corresponds to a given value and ``U`` to the applied voltage. For now,
the quantitity ``\\psi_S`` needs to be specified in the main file.
For the charge carriers we assume the following

``f[n_\\alpha]  =  z_\\alpha q v_\\alpha (n_\\alpha - n_{\\alpha, 0})``,

where ``v_{\\alpha}`` can be treated as a surface recombination mechanism and is given. The parameter
``n_{\\alpha, 0}`` is a given value, calculated by the statistical relation, when assuming 
no electrical field and a quasi Fermi level equal to the metal work function ``\\phi``, i.e.
    
``n_{\\alpha, 0}= z_\\alpha/ U_T (E_\\alpha - \\phi) / q. ``

"""
function breactionSchottky!(f, u, bnode, data)

    params        = data.params
    paramsnodal   = data.paramsnodal
    # DA: Still need to be well tested!!

    ipsi          = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1        # final index for electrostatic potential

    # bnode.coord
    # if bnode.region == 1
    #     f[ipsi] = -u[ipsi] 
    # elseif bnode.region == 2
    #     f[ipsi] =  + u[ipsi] - data.λ1 *((data.bFermiLevel[2] - data.bFermiLevel[1])/q ) - data.contactVoltage[bnode.region]  
    # end

    for icc = 1:params.numberOfCarriers-1
       
        if bnode.region == 1
            phi = params.bFermiLevel[1]
        elseif bnode.region == 2
            phi = params.bFermiLevel[2]
        end
        if bnode.region == 1 || bnode.region == 2
            E      = params.bBandEdgeEnergy[icc, bnode.region] + paramsnodal.bandEdgeEnergy[icc, bnode.index]
            etaFix = params.chargeNumbers[icc] / params.UT * (  (-phi + E ) / q  )
            eta    = params.chargeNumbers[icc] / params.UT * (  (u[icc]  - u[ipsi]) + E / q )

            f[icc] =  data.λ1 * params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * (  (params.bDensityOfStates[icc, bnode.region] + paramsnodal.densityOfStates[icc, bnode.index])  * (data.F[icc](eta) - data.F[icc](etaFix)  ))
        end

    end

    return f

end


##########################################################
##########################################################