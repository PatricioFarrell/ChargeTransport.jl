##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for interior nodes.
"""
function etaFunction(u, node::VoronoiFVM.Node, data, icc, ipsi)
     # params      = data.params
    # paramsnodal = data.paramsnodal

    E  = data.params.bandEdgeEnergy[icc, node.region] + data.paramsnodal.bandEdgeEnergy[icc, node.index]
    
     return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[ipsi]) + E / q )
   # return etaFunction(u[ipsi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc]) 
end

"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for boundary nodes.
"""
function etaFunction(u, bnode::VoronoiFVM.BNode, data, icc, ipsi) # bnode.index refers to index in overall mesh
    # params      = data.params
    # paramsnodal = data.paramsnodal

    E  = data.params.bBandEdgeEnergy[icc, bnode.region] + data.paramsnodal.bandEdgeEnergy[icc, bnode.index]
    
    #return params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
    return etaFunction(u[ipsi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc]) 
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for edges.
"""

function etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge; nodeside=1 )

    # params      = data.params
    # paramsnodal = data.paramsnodal
    
    E  = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, nodeEdge]

    # params.chargeNumbers[icc] / params.UT * ( (u[icc] - u[ipsi]) + E / q )
    return etaFunction(u[ipsi,nodeside], u[icc,nodeside], data.params.UT, E, data.params.chargeNumbers[icc]) 
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
function etaFunction(psi, phi, UT, E, z)
    z ./ UT .* ( (phi - psi) .+ E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for floats.
"""
function etaFunction(u, data, node, region, icc, ipsi, in_region::Bool)

    # params      = data.params
    # paramsnodal = data.paramsnodal

    if in_region == true
        E  = data.params.bandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    elseif in_region == false
        E  = data.params.bBandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    end

    return etaFunction(u[ipsi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc]) 
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
    ipsi  = data.indexPsi  # final index for electrostatic potential
 

    for icc ∈ data.chargeCarrierList # quantities or integer indices
 
        eta     = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )                    # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.bDensityOfStates[icc, bnode.region] + paramsnodal.densityOfStates[icc, bnode.index]) * data.F[icc](eta)  # add charge carrier
 
        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0
 
    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]

    f[ipsi] = - data.λ1 * 1 / tiny_penalty_value *  q * f[ipsi]

    return f

end


"""
$(TYPEDSIGNATURES)

This breaction! function is chosen when no interface model is chosen.

"""
breaction!(f, u, bnode, data, ::Type{interface_model_none}) = emptyFunction()


"""
$(TYPEDSIGNATURES)

breaction term for surface recombination.
"""

function breaction!(f, u, bnode, data, ::Type{interface_model_surface_recombination})
    if data.calculation_type == inEquilibrium
        return
    end

    # indices (∈ IN ) of electron and hole quasi Fermi potentials specified by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin # integer index of φ_n
    iphip       = data.bulk_recombination.iphip # integer index of φ_p

    ipsi        = data.indexPsi

    params      = data.params
    paramsnodal = data.paramsnodal
                                
    etan = params.chargeNumbers[iphin] / params.UT * ( (u[iphin] - u[ipsi]) + params.bBandEdgeEnergy[iphin, bnode.region] / q )
    etap = params.chargeNumbers[iphip] / params.UT * ( (u[iphip] - u[ipsi]) + params.bBandEdgeEnergy[iphip, bnode.region] / q ) 

    n    = ((params.bDensityOfStates[iphin, bnode.region] + paramsnodal.densityOfStates[iphin, bnode.index]) * FermiDiracZero(etan))
    p    = ((params.bDensityOfStates[iphip, bnode.region] + paramsnodal.densityOfStates[iphip, bnode.index]) * FermiDiracZero(etap))

    exponentialTerm = exp((q * u[iphin] - q  * u[iphip] ) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    kernelSRH = 1.0 / (  1.0/params.recombinationSRHvelocity[iphip, bnode.region] * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + 1.0/params.recombinationSRHvelocity[iphin, bnode.region] * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region] ) )
   
    for icc ∈ [iphin, iphip]
        f[icc] =  q * params.chargeNumbers[icc] * kernelSRH *  excessDensTerm
    end


end



"""
$(TYPEDSIGNATURES)

breaction term for case where qF are discontinuous.
"""

function breaction!(f, u, bnode, data, ::Type{interface_model_discont_qF})

    if data.calculation_type == inEquilibrium

        return emptyFunction()

    end

    ipsi = data.indexPsi

    #indices (∈ IN ) of electron and hole quasi Fermi potentials specified by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin # integer index of φ_n
    iphip       = data.bulk_recombination.iphip # integer index of φ_p

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin] # = Quantity or integer
    iphip       = data.chargeCarrierList[iphip] # = Quantity or integer

    params      = data.params
    paramsnodal = data.paramsnodal

    recombinationVelocity    = [1.0e1 1.0e7;
                                    1.0e5 1.0e1]
    
    ######### left values  ##########
    etan1 = params.chargeNumbers[iphin] / params.UT * ( (u[iphin, 1] - u[ipsi]) + params.bandEdgeEnergy[iphin, bnode.cellregions[1]] / q ) # left
    etap1 = params.chargeNumbers[iphip] / params.UT * ( (u[iphip, 1] - u[ipsi]) + params.bandEdgeEnergy[iphip, bnode.cellregions[1]] / q ) # left

    n1    = ((params.densityOfStates[iphin, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphin, bnode.index])* data.F[iphin](etan1))
    p1    = ((params.densityOfStates[iphip, bnode.cellregions[1]] + paramsnodal.densityOfStates[iphip, bnode.index])* data.F[iphip](etap1))

    exponentialTerm1 = exp((q * u[iphin, 1] - q  * u[iphip, 1] ) / (kB * params.temperature))
    excessDensTerm1  = n1 * p1 * (1.0 - exponentialTerm1)

    kernelSRH1 = 1.0 / (  1.0/recombinationVelocity[iphip, bnode.region-2] * (n1 + params.recombinationSRHTrapDensity[iphin, bnode.cellregions[1]]) + 1.0/recombinationVelocity[iphin, bnode.region-2] * (p1 + params.recombinationSRHTrapDensity[iphip, bnode.cellregions[1]] ) )


    ######### right values  ##########
    etan2 = params.chargeNumbers[iphin] / params.UT * ( (u[iphin, 2] - u[ipsi]) + params.bandEdgeEnergy[iphin, bnode.cellregions[2]] / q ) # right
    etap2 = params.chargeNumbers[iphip] / params.UT * ( (u[iphip, 2] - u[ipsi]) + params.bandEdgeEnergy[iphip, bnode.cellregions[2]] / q ) # right

    n2    = ((params.densityOfStates[iphin, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphin, bnode.index])* data.F[iphin](etan2))
    p2    = ((params.densityOfStates[iphip, bnode.cellregions[2]] + paramsnodal.densityOfStates[iphip, bnode.index])* data.F[iphip](etap2))

    exponentialTerm2 = exp((q * u[iphin, 2] - q  * u[iphip, 2] ) / (kB * params.temperature))
    excessDensTerm2   = n2 * p2 * (1.0 - exponentialTerm2)

    kernelSRH2 = 1.0 / (  1.0/recombinationVelocity[iphip, bnode.region-2] * (n2 + params.recombinationSRHTrapDensity[iphin, bnode.cellregions[2]])+ 1.0/recombinationVelocity[iphin, bnode.region-2] * (p2 + params.recombinationSRHTrapDensity[iphip, bnode.cellregions[2]]))

    for icc ∈ [iphin, iphip] # equations for qF potentials 
        react1     =  q *  params.chargeNumbers[icc] *  kernelSRH1 *  excessDensTerm1 
        react2     =  q *  params.chargeNumbers[icc] *  kernelSRH2 *  excessDensTerm2 

        f[icc, 1] = react1
        f[icc, 2] = react2 # same plot with minus sign here

    end

    # for icc ∈ [iphin, iphip] # equations for qF potentials 
            
    #     #to see continuity
    #     d         = [1.0e7 1.0e7;
    #                1.0e7 1.0e7]

    #     # to see discontinuity
    #     # d         = [1.0e1 1.0e3;
    #     #                 1.0e7 1.0e1]
    #     react     =d[icc, bnode.region-2] *   (u[icc, 1] - u[icc, 2])

    #     f[icc, 1] =   react
    #     f[icc, 2] = - react
    # end

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
No bstorage! is used, when assuming discontinuous qF.

"""
bstorage!(f, u, bnode, data, ::Type{interface_model_discont_qF}) = emptyFunction()

"""
$(TYPEDSIGNATURES)
No bstorage! is used, if an ohmic contact model is chosen.

"""
bstorage!(f, u, bnode, data, ::Type{ohmic_contact}) = emptyFunction()

"""
$(TYPEDSIGNATURES)
No bstorage! is used, if an schottky contact model is chosen.

"""
bstorage!(f, u, bnode, data, ::Type{schottky_contact}) = emptyFunction()



"""
$(TYPEDSIGNATURES)
Time-dependent part in case of present ion charges at inner interfaces (left).

"""
function bstorage!(f, u, bnode, data, ::Type{interface_model_ion_charge_left})
# DA: I guess, here is no sign included so we do not even need to distinguish between left and right.
# But need to fix use of interface charges (done by indices)!

    params            = data.params
    paramsnodal       = data.paramsnodal

    iphia             = data.chargeCarrierList[3]
    iphiaj1, iphiaj2  = 4:5
    ipsi              = data.indexPsi      # final index for electrostatic potential

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
Master bflux! function. This is the function which enters VoronoiFVM and hands over
for each boundary the flux within the boundary.

"""
bflux!(f, u, bedge, data) = bflux!(f, u, bedge, data, data.calculation_type)


"""
In case of equilibrium, the bflux shall not enter.
"""
bflux!(f, u, bedge, data, ::Type{inEquilibrium}) = emptyFunction()



"""
Out of equilibrium, we need to additionally check for grid dimension.
"""
bflux!(f, u, bedge, data, ::Type{outOfEquilibrium}) = bflux!(f, u, bedge, data, data.grid_dimension)


"""
In case of one dimensional grid, no bflux entering.
"""
bflux!(f, u, bedge, data, ::Type{OneD_grid}) = emptyFunction()


"""
Out of equilibrium and for dimension = 2, the bflux shall only enter, when we have inner interfaces defined.
"""
bflux!(f, u, bedge, data, ::Type{TwoD_grid}) = bflux!(f, u, bedge, data, data.boundary_type[bedge.region]) #emptyFunction()#


"""
For outer boundaries.
"""
bflux!(f, u, bedge, data, ::Type{ohmic_contact})    = emptyFunction()
bflux!(f, u, bedge, data, ::Type{schottky_contact}) = emptyFunction()

"""
In this specific case then, we can use the precise flux approximation scheme.

"""
# DA: does not work with Type{interface_model} only ????
bflux!(f, u, bedge, data, ::Type{interface_model_none}) = bflux!(f, u, bedge, data, data.flux_approximation)
bflux!(f, u, bedge, data, ::Type{interface_model_surface_recombination}) = bflux!(f, u, bedge, data, data.flux_approximation)

"""
$(TYPEDSIGNATURES)

The excess chemical potential flux discretization scheme for inner boundaries.

"""
function bflux!(f, u, bedge, data, ::Type{excessChemicalPotential})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    
    ipsi        =   data.indexPsi
    nodel       =   bedge.node[2]
    nodek       =   bedge.node[1]
    ireg        =   bedge.region
    
    # ############################################################
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
  
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList[1:2]

        j0                 = params.chargeNumbers[icc] * q * params.bMobility[icc, ireg] * params.UT * params.bDensityOfStates[icc, ireg]

        # need to add this to the other etaFunctions
        Ek                 = params.bBandEdgeEnergy[icc, bedge.region] + paramsnodal.bandEdgeEnergy[icc, nodek]
        El                 = params.bBandEdgeEnergy[icc, bedge.region] + paramsnodal.bandEdgeEnergy[icc, nodel]
        etak               = etaFunction(u[ipsi, 1], u[icc, 1], data.params.UT, Ek, params.chargeNumbers[icc])
        etal               = etaFunction(u[ipsi, 2], u[icc, 2], data.params.UT, El, params.chargeNumbers[icc])

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q                  = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )

    end

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
    ireg        = node.region
    inode       = node.index

    ipsi        = data.indexPsi   # final index for electrostatic potential

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, ireg] )   # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, inode]) * data.F[icc](eta)   # add charge carrier

    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[inode]

    f[ipsi] = - data.λ1 * q * f[ipsi]

    ############################################################
    ####          simple zero reaction for all icc          ####
    ############################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}

        f[icc] = u[icc] - 0.0

    end

    return f

end

recombination_kernel(data, ireg, iphin, iphip, n, p, ::Type{bulk_recomb_model_none}) = 0.0


function recombination_kernel(data, ireg, iphin, iphip,  n, p, ::Type{bulk_recomb_model_radiative})

    params = data.params

    return params.recombinationRadiative[ireg]

end


function recombination_kernel(data, ireg, iphin, iphip, n, p,::Type{bulk_recomb_model_trap_assisted})

    params = data.params

    kernelSRH = 1.0 / (  params.recombinationSRHLifetime[iphip, ireg] * (n + params.recombinationSRHTrapDensity[iphin, ireg]) + params.recombinationSRHLifetime[iphin, ireg] * (p + params.recombinationSRHTrapDensity[iphip, ireg]) )

    return kernelSRH

end


function recombination_kernel(data, ireg, iphin, iphip, n, p, ::Type{bulk_recomb_model_full})

    params = data.params

    # radiative recombination
    kernelRadiative = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recomb_model_radiative)

    # SRH recombination
    kernelSRH       = recombination_kernel(data, ireg, iphin, iphip, n, p, bulk_recomb_model_trap_assisted)

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
    ireg        = node.region
    inode       = node.index

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.indexPsi                  # final index for electrostatic potential
    
    ############################################################
    ####          simple zero reaction for all icc          ####
    ############################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}

        f[icc]  = u[icc] - 0.0 # set for all charge carriers (electric and possible present ionic) zero conditions

    end

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    f[ipsi]     = f[ipsi] - paramsnodal.doping[inode]

    f[ipsi]     = - q * data.λ1 * f[ipsi]

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    n               = (params.densityOfStates[iphin, ireg] + paramsnodal.densityOfStates[iphin, inode]) * data.F[iphin](etaFunction(u, node, data, iphin, ipsi))
    p               = (params.densityOfStates[iphip, ireg] + paramsnodal.densityOfStates[iphip, inode]) * data.F[iphip](etaFunction(u, node, data, iphip, ipsi))

    exponentialTerm = exp((q *u[iphin] - q * u[iphip]) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    for icc ∈ [iphin, iphip] 

        # gives you the recombination kernel based on choice of user
        kernel      = recombination_kernel(data, ireg, iphin, iphip, n, p, data.bulk_recombination.bulk_recomb_model)
                
        f[icc]      = q * params.chargeNumbers[icc] *  kernel *  excessDensTerm  - q * params.chargeNumbers[icc] * generation(data, ireg,  node.coord[node.index], data.generation_model)
    end

    return f

end

"""
$(TYPEDSIGNATURES)

Like ``reaction!(f, u, node, data, ::Type{outOfEquilibrium})`` but including a trap density for transient simulations.

"""
function reaction!(f, u, node, data, ::Type{outOfEquilibrium_trap})

    params      = data.params
    paramsnodal = data.paramsnodal
    ireg        = node.region
    inode       = node.index

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used 
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    itrap       = data.chargeCarrierList[3]
    ipsi        = data.indexPsi                  

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    f[ipsi]     = f[ipsi] - paramsnodal.doping[inode]
    f[ipsi]     = - q * f[ipsi]

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    Nc = params.densityOfStates[iphin, ireg] + paramsnodal.densityOfStates[iphin, inode]
    Nv = params.densityOfStates[iphip, ireg] + paramsnodal.densityOfStates[iphip, inode]
    Nt = params.densityOfStates[itrap, ireg]
    n  = Nc * data.F[iphin](etaFunction(u, node, data, iphin, ipsi))
    p  = Nv * data.F[iphip](etaFunction(u, node, data, iphip, ipsi))
    t  = Nt * data.F[itrap](etaFunction(u, node, data, itrap, ipsi))

    taun                  = params.recombinationSRHLifetime[iphin, ireg]
    n0                    = params.recombinationSRHTrapDensity[iphin, ireg]
    taup                  = params.recombinationSRHLifetime[iphip, ireg]
    p0                    = params.recombinationSRHTrapDensity[iphip, ireg]

    # Rn, Rp agree up to sign with *On the Shockley-Read-Hall Model: Generation-Recombination 
    # in Semiconductors* in SIAM Journal on Applied Mathematics, Vol. 67, No. 4 (2007), pp. 1183-1201.
    # The sign is chosen according to *Supporting Information: Consistent Device Simulation Model 
    # Describing Perovskite Solar Cells in Steady-State, Transient and Frequency Domain* in ACS (2018)
    if params.chargeNumbers[itrap] == -1
        Rn =  1 / taun * (n * (1-t/Nt) - n0 * t/Nt)
        Rp =  1 / taup * (p * t/Nt     - p0 * (1-t/Nt))
    elseif params.chargeNumbers[itrap] == 1
        Rn =  1 / taun * (n * t/Nt     - n0 * (1-t/Nt))
        Rp =  1 / taup * (p * (1-t/Nt) - p0 * t/Nt)
    end

    #@show generation(data, ireg,  node.coord[node.index], data.generation_model)
    f[iphin] = q * params.chargeNumbers[iphin] * (Rn - generation(data, ireg,  node.coord[node.index], data.generation_model))
    f[iphip] = q * params.chargeNumbers[iphip] * (Rp - generation(data, ireg,  node.coord[node.index], data.generation_model))
    f[itrap] = q * params.chargeNumbers[itrap] * (Rp-Rn)

    return f    
end

"""
$(TYPEDSIGNATURES)

Like ``reaction!(f, u, node, data, ::Type{outOfEquilibrium_trap})`` but for stationary simulations.

"""
function reaction!(f, u, node, data, ::Type{outOfEquilibrium_trap_stationary})

    params      = data.params
    paramsnodal = data.paramsnodal
    ireg        = node.region
    inode       = node.index

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    itrap       = 3
    ipsi        = data.indexPsi                  # final index for electrostatic potential
    
    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################

    for icc ∈ data.chargeCarrierList 
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    f[ipsi]     = f[ipsi] - paramsnodal.doping[inode]
    f[ipsi]     = - q * f[ipsi]

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    n  = (params.densityOfStates[iphin, ireg] + paramsnodal.densityOfStates[iphin, inode]) * data.F[iphin](etaFunction(u, node, data, iphin, ipsi))
    p  = (params.densityOfStates[iphip, ireg] + paramsnodal.densityOfStates[iphip, inode]) * data.F[iphip](etaFunction(u, node, data, iphip, ipsi))

    exponentialTerm = exp((q *u[iphin] - q * u[iphip]) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    for icc ∈ [iphin, iphip] 

        # gives you the recombination kernel based on choice of user
        kernel      = recombination_kernel(data, ireg, iphin, iphip, n, p, data.bulk_recombination.bulk_recomb_model)
                
        f[icc]      = q * params.chargeNumbers[icc] *  kernel *  excessDensTerm  - q * params.chargeNumbers[icc] * generation(data, ireg,  node.coord[node.index], data.generation_model)
    end

    # switch off traps by inserting unit matrix into Jacobian
    f[itrap]  = u[itrap] - 0.0

    return f

end

"""
$(TYPEDSIGNATURES)

Like ``reaction!(f, u, node, data, ::Type{outOfEquilibrium_trap})`` but for stationary simulations and only two species. The trap density is explicitly added to the right-hand side.

"""
function reaction!(f, u, node, data, ::Type{outOfEquilibrium_trap_stationary_2_species})

    params      = data.params
    paramsnodal = data.paramsnodal
    ireg        = node.region
    inode       = node.index

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.indexPsi                  # final index for electrostatic potential
    
    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################

    for icc ∈ data.chargeCarrierList 
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * (params.densityOfStates[icc, node.region] + paramsnodal.densityOfStates[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    n  = (params.densityOfStates[iphin, ireg] + paramsnodal.densityOfStates[iphin, inode]) * data.F[iphin](etaFunction(u, node, data, iphin, ipsi))
    p  = (params.densityOfStates[iphip, ireg] + paramsnodal.densityOfStates[iphip, inode]) * data.F[iphip](etaFunction(u, node, data, iphip, ipsi))
    n0 = params.recombinationSRHTrapDensity[iphin, ireg]
    p0 = params.recombinationSRHTrapDensity[iphip, ireg]
    taun = params.recombinationSRHLifetime[iphin, ireg]
    taup = params.recombinationSRHLifetime[iphip, ireg]
    z = 1

    if ireg == 3
        Nt = 5e14                / (cm^3) 
    else
        Nt = 5e14                / (cm^3) 
    end

    # add equilibrium trap density
    f[ipsi] = f[ipsi] + z * Nt * (taun*p0 + taup*n) / (taun*(p0+p) + taup*(n0+n))


    f[ipsi]     = f[ipsi] - paramsnodal.doping[inode]
    f[ipsi]     = - q * f[ipsi]

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    exponentialTerm = exp((q *u[iphin] - q * u[iphip]) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    for icc ∈ [iphin, iphip] 

        # gives you the recombination kernel based on choice of user
        kernel      = recombination_kernel(data, ireg, iphin, iphip, n, p, data.bulk_recombination.bulk_recomb_model)
                
        f[icc]      = q * params.chargeNumbers[icc] *  kernel *  excessDensTerm  - q * params.chargeNumbers[icc] * generation(data, ireg,  node.coord[node.index], data.generation_model)
    end

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

    ipsi        = data.indexPsi
    
    for icc ∈ data.chargeCarrierList

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
    paramsnodal = data.paramsnodal
    
    ipsi        = data.indexPsi # integer index or Quantity
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - (params.dielectricConstant[ireg] + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * ε0 * dpsi
    
end

flux!(f, u, edge, data, ::Type{outOfEquilibrium}) = flux!(f, u, edge, data, data.flux_approximation)

function flux!(f, u, edge, data, ::Type{outOfEquilibrium_trap}) 
    flux!(f, u, edge, data, outOfEquilibrium) 
end

function flux!(f, u, edge, data, ::Type{outOfEquilibrium_trap_stationary}) 
    flux!(f, u, edge, data, outOfEquilibrium) 
end

function flux!(f, u, edge, data, ::Type{outOfEquilibrium_trap_stationary_2_species}) 
    flux!(f, u, edge, data, outOfEquilibrium) 
end


flux!(f, u, edge, data, ::Type{flux_approximation}) = emptyFunction()

"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{ScharfetterGummel})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.indexPsi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 =  params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2)

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

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.indexPsi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 =  params.chargeNumbers[icc] * q * params.UT

        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2) 

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobility           = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2
        densityOfStatesl   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, nodel])
        densityOfStatesk   = (params.densityOfStates[icc, ireg] + paramsnodal.densityOfStates[icc, nodek])
        
        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT ) 
        else
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc, nodek])) ) 
        end

        f[icc]             = - j0  * mobility * ( bm  * densityOfStatesl * data.F[icc](etal) - bp *  densityOfStatesk * data.F[icc](etak) )

    end

    return f

end


"""
$(TYPEDSIGNATURES)

The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
not for space-dependent effective DOS.

"""
function flux!(f, u, edge, data, ::Type{excessChemicalPotential})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.indexPsi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2) 

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q                  = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
        # @show value(data.F[icc](etal))
        # @show value(data.F[icc](etak))
        # @show value(Q)
        # @show value((dpsi - bandEdgeDifference/q) /params.UT)
        # @show value(u[3])
    end


    #     @show value(f[1])
    # @show value(f[2])
    # @show value(f[3])


    
    # println("---")

    return f

end


"""
$(TYPEDSIGNATURES)

The excess chemical potential flux scheme for 
possible space-dependent DOS and band-edge energies. For these parameters
the discretization scheme is modified. [insert continuous flux etc ...]

"""
function flux!(f, u, edge, data, ::Type{excessChemicalPotential_Graded})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.indexPsi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.UT
 
        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2) 

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

    params      =   data.params
    paramsnodal =   data.paramsnodal

    tolReg      =   1.0e-13
    
    ipsi        =   data.indexPsi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2)

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

    params        =   data.params
    paramsnodal   =   data.paramsnodal

    max_iteration =   200          # for Newton solver
    it            =   0            # number of iterations (newton)
    damp          =   0.1          # damping factor
    
    ipsi          =   data.indexPsi
    nodel         =   edge.node[2]
    nodek         =   edge.node[1]
    ireg          =   edge.region
    
    dpsi          =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]       =  - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                  = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT * params.densityOfStates[icc, ireg]

        # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge, nodeside)
        etak               = etaFunction(u, edge, data, icc, ipsi, nodek, nodeside = 1) 
        etal               = etaFunction(u, edge, data, icc, ipsi, nodel, nodeside = 2)

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
function breaction!(f, u, bnode, data,  ::Type{schottky_contact})

    params        = data.params
    paramsnodal   = data.paramsnodal

    ipsi          = params.numberOfCarriers + params.numberOfInterfaceCarriers + 1        # final index for electrostatic potential

    for icc = 1:params.numberOfCarriers
       
        E      = params.bBandEdgeEnergy[icc, bnode.region] + paramsnodal.bandEdgeEnergy[icc, bnode.index]
        etaFix = params.chargeNumbers[icc] / params.UT * (  (- params.bFermiLevel[bnode.region] + E ) / q  )
        eta    = params.chargeNumbers[icc] / params.UT * (  (u[icc]  - u[ipsi]) + E / q )

        f[icc] =  data.λ1 * params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * (  (params.bDensityOfStates[icc, bnode.region] + paramsnodal.densityOfStates[icc, bnode.index])  * (data.F[icc](eta) - data.F[icc](etaFix)  ))

    end

    return f

end


##########################################################
##########################################################