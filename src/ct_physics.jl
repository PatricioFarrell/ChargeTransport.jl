##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the effective DOS for interior nodes (analogously for boundary nodes and edges).
"""
function get_DOS!(icc, node::VoronoiFVM.Node, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, node.region] + data.paramsnodal.densityOfStates[icc, node.index]

end

# Defining locally the effective DOS for boundary nodes.
function get_DOS!(icc, bnode::VoronoiFVM.BNode, data)

    data.tempDOS1[icc] = data.params.bDensityOfStates[icc, bnode.region] + data.paramsnodal.densityOfStates[icc, bnode.index]

end

# Defining locally the effective DOS for edges.
function get_DOS!(icc, edge::VoronoiFVM.Edge, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[1]]
    data.tempDOS2[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[2]]

    return data.tempDOS1, data.tempDOS2
end

"""
$(TYPEDSIGNATURES)

Defining locally the band-edge energy for interior nodes (analougesly for boundary nodes and edges).
"""
function get_BEE!(icc, node::VoronoiFVM.Node, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, node.region] + data.paramsnodal.bandEdgeEnergy[icc, node.index]

end

# Defining locally the band-edge energy for boundary nodes.
function get_BEE!(icc, bnode::VoronoiFVM.BNode, data)

    data.tempBEE1[icc] = data.params.bBandEdgeEnergy[icc, bnode.region] + data.paramsnodal.bandEdgeEnergy[icc, bnode.index]

end

# Defining locally the band-edge energy for edges.
function get_BEE!(icc, edge::VoronoiFVM.Edge, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[1]]
    data.tempBEE2[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[2]]
    
end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the distribution function for interior nodes.
"""
function etaFunction(u, node::VoronoiFVM.Node, data, icc, ipsi)

    get_BEE!(icc, node::VoronoiFVM.Node, data)

    E  = data.tempBEE1[icc]
    
    return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[ipsi]) + E / q )

end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for floats.
"""
function etaFunction(u, data, node, region, icc, ipsi, in_region::Bool)

    if in_region == true
        E  = data.params.bandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    elseif in_region == false
        E  = data.params.bBandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    end

    return etaFunction(u[ipsi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc]) 
end

"""
$(TYPEDSIGNATURES)

The argument of the distribution function for boundary nodes.
"""
function etaFunction(u, bnode::VoronoiFVM.BNode, data, icc, ipsi) # bnode.index refers to index in overall mesh

    get_BEE!(icc, bnode::VoronoiFVM.BNode, data)

    E  = data.tempBEE1[icc]
    
    return etaFunction(u[ipsi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc]) 
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function for edges.
"""

function etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

    get_BEE!(icc, edge::VoronoiFVM.Edge, data)

    E1 = data.tempBEE1[icc];  E2 = data.tempBEE2[icc]

    return etaFunction(u[ipsi, 1], u[icc, 1], data.params.UT, E1, data.params.chargeNumbers[icc]), etaFunction(u[ipsi, 2], u[icc, 2], data.params.UT, E2, data.params.chargeNumbers[icc]) 
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
    ipsi        = data.index_psi  # final index for electrostatic potential
 
    for icc ∈ data.chargeCarrierList # quantities or integer indices
 
        get_DOS!(icc, bnode, data)
        Ni      = data.tempDOS1[icc]
        eta     = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )   # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * Ni * data.F[icc](eta)                   # add charge carrier
 
        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0
 
    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]

    f[ipsi] = - data.λ1 * 1 / tiny_penalty_value *  q * f[ipsi]

    return f

end

"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions in a first attempt. For the electrostatic potential we assume 

``\\psi = \\phi_S/q + U, ``

where  ``\\phi_S`` corresponds to a given value (Schottky barrier) and ``U`` to the applied voltage. For now,
the quantitity ``\\phi_S`` needs to be specified in the main file.
For the charge carriers we assume the following

``f[n_\\alpha]  =  z_\\alpha q v_\\alpha (n_\\alpha - n_{\\alpha, 0})``,

where ``v_{\\alpha}`` can be treated as a surface recombination mechanism and is given. The parameter
``n_{\\alpha, 0}`` is a given value, calculated by the statistical relation, when assuming 
no electrical field and a quasi Fermi level equal to the Schottky barrier ``\\phi_S``, i.e.
    
``n_{\\alpha, 0}= z_\\alpha/ U_T (E_\\alpha - \\phi_S) / q. ``

"""
function breaction!(f, u, bnode, data,  ::Type{schottky_contact})

    params        = data.params
    paramsnodal   = data.paramsnodal

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used 
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.index_psi               

    for icc in [iphin,iphip] 

        get_DOS!(icc, bnode, data);  get_BEE!(icc, bnode, data)
        Ni     = data.tempDOS1[icc]
        Ei     = data.tempBEE1[icc]
        etaFix = params.chargeNumbers[icc] / params.UT * (  (- params.SchottkyBarrier[bnode.region] + Ei ) / q  )
        eta    = params.chargeNumbers[icc] / params.UT * (  (u[icc]  - u[ipsi]) + Ei / q )

        f[icc] =  data.λ1 * params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * (  Ni  * (data.F[icc](eta) - data.F[icc](etaFix)  ))

    end

    return f

end


# This breaction! function is chosen when no interface model is chosen.
breaction!(f, u, bnode, data, ::Type{interface_model_none}) = emptyFunction()


# breaction term for surface recombination.
breaction!(f, u, bnode, data, ::Type{interface_model_surface_recombination_and_tangential_flux}) = breaction!(f, u, bnode, data, interface_model_surface_recombination)


function breaction!(f, u, bnode, data, ::Type{interface_model_surface_recombination})
    if data.calculation_type == inEquilibrium
        return
    end

    # indices (∈ IN ) of electron and hole quasi Fermi potentials specified by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin # integer index of φ_n
    iphip       = data.bulk_recombination.iphip # integer index of φ_p

    ipsi        = data.index_psi

    params      = data.params
    paramsnodal = data.paramsnodal
    
    get_DOS!(iphin, bnode, data); get_DOS!(iphip, bnode, data)
    Nc   = data.tempDOS1[iphin]
    Nv   = data.tempDOS1[iphip]
    etan = etaFunction(u, bnode, data, iphin, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
    etap = etaFunction(u, bnode, data, iphip, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)

    n    = Nc * data.F[iphin](etan)
    p    = Nv * data.F[iphip](etap)

    exponentialTerm = exp((q * u[iphin] - q  * u[iphip] ) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    kernelSRH = 1.0 / (  1.0/params.recombinationSRHvelocity[iphip, bnode.region] * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + 1.0/params.recombinationSRHvelocity[iphin, bnode.region] * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region] ) )
   
    for icc ∈ [iphin, iphip]
        f[icc] =  q * params.chargeNumbers[icc] * kernelSRH *  excessDensTerm
    end


end

# breaction term for case where qF are discontinuous.
function breaction!(f, u, bnode, data, ::Type{interface_model_discont_qF})

    if data.calculation_type == inEquilibrium

        return emptyFunction()

    end

    ipsi = data.index_psi

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
bstorage!(f, u, bnode, data) = bstorage!(f, u, bnode, data, data.model_type)



bstorage!(f, u, bnode, data, ::Type{model_stationary})  = emptyFunction()


bstorage!(f, u, bnode, data, ::Type{model_transient}) = bstorage!(f, u, bnode, data, data.boundary_type[bnode.region])


bstorage!(f, u, bnode, data, ::Type{interface_model_none}) = emptyFunction()


# No bstorage! is used, when assuming discontinuous qF.
bstorage!(f, u, bnode, data, ::Type{interface_model_discont_qF}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{interface_model_surface_recombination}) = emptyFunction()

# No bstorage! is used, if an ohmic and schottky contact model is chosen.
bstorage!(f, u, bnode, data, ::Type{ohmic_contact}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{schottky_contact}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{interface_model_surface_recombination_and_tangential_flux}) = bstorage!(f, u, bnode, data, interface_model_tangential_flux)

function bstorage!(f, u, bnode, data, ::Type{interface_model_tangential_flux})

    params      = data.params
    paramsnodal = data.paramsnodal 

    #indices (∈ IN ) of electron and hole quasi Fermi potentials specified by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin # integer index of φ_n
    iphip       = data.bulk_recombination.iphip # integer index of φ_p

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin] # = Quantity or integer
    iphip       = data.chargeCarrierList[iphip] # = Quantity or integer

    ipsi        = data.index_psi
    
    for icc ∈ [iphin, iphip]

        get_DOS!(icc, bnode, data)
        Ni     = data.tempDOS[icc]
        eta    = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = q * params.chargeNumbers[icc] * Ni * data.F[icc](eta)

    end

    f[ipsi] =  0.0

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


# In case of equilibrium, the bflux shall not enter.
bflux!(f, u, bedge, data, ::Type{inEquilibrium})                         = emptyFunction()


# Out of equilibrium, we need to additionally check for boundary type.
bflux!(f, u, bedge, data, ::Type{outOfEquilibrium})                      = bflux!(f, u, bedge, data, data.boundary_type[bedge.region])

bflux!(f, u, bedge, data, ::Type{interface_model_none})                  = emptyFunction()

bflux!(f, u, bedge, data, ::Type{ohmic_contact})                         = emptyFunction()
bflux!(f, u, bedge, data, ::Type{schottky_contact})                      = emptyFunction()
bflux!(f, u, bedge, data, ::Type{interface_model_surface_recombination}) = emptyFunction()



# Cases, where we have a tangential flux.
bflux!(f, u, bedge, data, ::Type{interface_model_surface_recombination_and_tangential_flux}) = bflux!(f, u, bedge, data, data.flux_approximation)

bflux!(f, u, bedge, data, ::Type{interface_model_tangential_flux}) = bflux!(f, u, bedge, data, data.flux_approximation)


# excess chemical potential flux discretization scheme for inner boundaries.
function bflux!(f, u, bedge, data, ::Type{excess_chemical_potential})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    
    nodel       =   bedge.node[2]
    nodek       =   bedge.node[1]
    ireg        =   bedge.region

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.index_psi                  # final index for electrostatic potential
    
    # ############################################################
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
  
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ [iphin, iphip]

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

    ipsi        = data.index_psi  # final index for electrostatic potential

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        get_DOS!(icc, node, data)
        Ni      = data.tempDOS1[icc]
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.doping[icc, ireg] )   # subtract doping
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * Ni * data.F[icc](eta)   # add charge carrier

    end
    f[ipsi] = f[ipsi] - paramsnodal.doping[inode]

    f[ipsi] = - data.λ1 * q * f[ipsi]

    ############################################################
    ####            zero reaction term for all icc          ####
    ############################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}
        f[icc] = u[icc]
    end

    return f

end


"""
$(TYPEDSIGNATURES)
SRH kernel for case of non-existing rate.

"""
function kernelSRH(data, ireg, iphin, iphip, n, p, ::Type{model_SRH_off})

    return 0.0

end

"""
$(TYPEDSIGNATURES)
SRH kernel for case of using stationary formula, i.e. case where no present traps are assumed.

"""
function kernelSRH(data, ireg, iphin, iphip, n, p, ::Type{model_SRH_stationary})

    return  1.0 / (  data.params.recombinationSRHLifetime[iphip, ireg] * (n + data.params.recombinationSRHTrapDensity[iphin, ireg]) + data.params.recombinationSRHLifetime[iphin, ireg] * (p + data.params.recombinationSRHTrapDensity[iphip, ireg]) )

end



"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of Poisson equation, i.e. which builds
the space charge density for outOfEquilibrium calculations.

"""
function addElectricPotential!(f, u, node, data, ipsi)

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}

        get_DOS!(icc, node, data)

        Ni      = data.tempDOS1[icc]
        eta     = etaFunction(u, node, data, icc, ipsi) 

        f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * ( data.params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * Ni * data.F[icc](eta)   # add charge carrier

    end

    f[ipsi]     = f[ipsi] - data.paramsnodal.doping[node.index]

    f[ipsi]     = - q * f[ipsi]

    return f

end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of charge carriers.

"""
function addChargeCarriers!(f, u, node, data, ipsi, iphin, iphip, n, p)

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    # dependent on user information concerncing recombination, RHS will be build
    addRecombinationProcess!(f, u, node, data, ipsi, iphin, iphip, n, p, data.bulk_recombination.bulk_recomb_SRH)
   
    return f

end

function addRecombinationProcess!(f, u, node, data, ipsi, iphin, iphip, n, p, ::Type{T}) where T<:model_SRH_without_traps

    ireg        = node.region

    exponentialTerm = exp((q *u[iphin] - q * u[iphip]) / (kB * data.params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    for icc ∈ [iphin, iphip] 

        # gives you the recombination kernel based on choice of user
        # by default parameters of Auger and radiative are 0.0. Hence, adding them here, has no influence since
        # we simply add by 0.0.
        kernel      = data.params.recombinationRadiative[ireg] + (data.params.recombinationAuger[iphin, ireg] * n + data.params.recombinationAuger[iphip, ireg] * p)
        kernel      = kernel + kernelSRH(data, ireg, iphin, iphip, n, p, data.bulk_recombination.bulk_recomb_SRH)
                
        f[icc]      = q * data.params.chargeNumbers[icc] *  kernel *  excessDensTerm  - q * data.params.chargeNumbers[icc] * generation(data, ireg,  node.coord[node.index], data.generation_model)
    end

    return f
end

function addRecombinationProcess!(f, u, node, data, ipsi, iphin, iphip, n, p, ::Type{model_SRH_traps_transient})

    params      = data.params
    ireg        = node.region

    # indices (∈ IN ) of traps used by user (they pass it through recombination)
    itrap       = data.enable_traps.traps

    # based on user index and regularity of solution quantities or integers are used 
    itrap       = data.chargeCarrierList[itrap]

    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################

    Nt    = params.densityOfStates[itrap, ireg]
    t     = Nt * data.F[itrap](etaFunction(u, node, data, itrap, ipsi))

    taun  = params.recombinationSRHLifetime[iphin, ireg]
    n0    = params.recombinationSRHTrapDensity[iphin, ireg]
    taup  = params.recombinationSRHLifetime[iphip, ireg]
    p0    = params.recombinationSRHTrapDensity[iphip, ireg]

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

    f[iphin] = q * params.chargeNumbers[iphin] * (Rn - generation(data, ireg,  node.coord[node.index], data.generation_model))
    f[iphip] = q * params.chargeNumbers[iphip] * (Rp - generation(data, ireg,  node.coord[node.index], data.generation_model))
    f[itrap] = q * params.chargeNumbers[itrap] * (Rp-Rn)


    return f
end



# Function which adds additional trap density to right-hand side of Poisson equation
# without modeling traps as own charge carrier.
# Note that this one may be deleted in future version.

addTrapDensity!(f, u, node, data, ipsi, iphin, iphip, n, p) = addTrapDensity!(f, u, node, data, ipsi, iphin, iphip, n, p, data.bulk_recombination.model_SRH_2species_trap)


function addTrapDensity!(f, u, node, data, ipsi, iphin, iphip, n, p, ::Type{T}) where T<:model_SRH
    return
end

function addTrapDensity!(f, u, node, data, ipsi, iphin, iphip, n, p, ::Type{model_SRH_2species_present_trap_dens})

    params = data.params
    ireg   = node.region

    n0     = params.recombinationSRHTrapDensity[iphin, ireg]
    p0     = params.recombinationSRHTrapDensity[iphip, ireg]
    taun   = params.recombinationSRHLifetime[iphin, ireg]
    taup   = params.recombinationSRHLifetime[iphip, ireg]
    z      = 1

    if ireg == 3
        Nt = 5e14                / (cm^3) 
    else
        Nt = 5e14                / (cm^3) 
    end

    # add equilibrium trap density
    f[ipsi] = f[ipsi] - q * z * Nt * (taun*p0 + taup*n) / (taun*(p0+p) + taup*(n0+n))

    return f
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

    # indices (∈ IN ) of electron and hole quasi Fermi potentials used by user (they pass it through recombination)
    iphin       = data.bulk_recombination.iphin
    iphip       = data.bulk_recombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin       = data.chargeCarrierList[iphin]
    iphip       = data.chargeCarrierList[iphip]
    ipsi        = data.index_psi                # final index for electrostatic potential
    
    ############################################################
    ####   set RHS to zero for all icc (stability purpose)  ####
    ############################################################
    for icc ∈ data.chargeCarrierList # chargeCarrierList[icc] ∈ {IN} ∪ {AbstractQuantity}

        f[icc]  = u[icc] # set for all charge carriers right hand-side to zero

    end

    get_DOS!(iphin, node, data); get_DOS!(iphip, node, data)

    Nc          = data.tempDOS1[iphin];      Nv          = data.tempDOS1[iphip]

    n           = Nc * data.F[iphin](etaFunction(u, node, data, iphin, ipsi))
    p           = Nv * data.F[iphip](etaFunction(u, node, data, iphip, ipsi))

    
    f = addElectricPotential!(f, u, node, data, ipsi)                  # RHS of Poisson
    f = addChargeCarriers!(f, u, node, data, ipsi, iphin, iphip, n, p) # RHS of Charge Carriers with special treatment of recombination
    f = addTrapDensity!(f, u, node, data, ipsi, iphin, iphip, n, p)    # if desired, add trap density to RHS of Poisson

    return f

end


"""
$(SIGNATURES)

Compute trap densities for a given trap energy.
[Currently, only done for the Boltzmann statistics and for region dependent parameters.]

"""
function trap_density!(icc, ireg, data, Et) 
    params      = data.params

    params.densityOfStates[icc, ireg] * exp( params.chargeNumbers[icc] * (params.bandEdgeEnergy[icc, ireg] - Et) / (kB * params.temperature))
end

# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations with a uniform generation rate.
function generation(data, ireg, node, ::Type{generation_uniform}) # only works in 1D till now; adjust node, when multidimensions

    params = data.params

    return data.λ2 * params.generationUniform[ireg]
end


# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations obeying the Beer-Lambert law.
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

    ipsi        = data.index_psi
    
    for icc ∈ data.chargeCarrierList

        get_DOS!(icc, node, data)

        Ni     = data.tempDOS1[icc]
        eta    = etaFunction(u, node, data, icc, ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = q * params.chargeNumbers[icc] * Ni * data.F[icc](eta)

    end

    f[ipsi] =  0.0

    return f
end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master flux functions which enters VoronoiFVM. Flux discretization scheme is chosen in two steps. First, we need
to see, if we are in or out of equilibrium. If, inEquilibrium, then
no flux is passed. If outOfEquilibrium, we choose the flux approximation
which the user chose.

"""
flux!(f, u, edge, data) = flux!(f, u, edge, data, data.calculation_type)

function flux!(f, u, edge, data, ::Type{inEquilibrium})

    params      = data.params
    paramsnodal = data.paramsnodal
    
    ipsi        = data.index_psi # integer index or Quantity
    ireg        = edge.region
    nodel       = edge.node[2]
    nodek       = edge.node[1]

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - (params.dielectricConstant[ireg] + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * ε0 * dpsi
    
end

flux!(f, u, edge, data, ::Type{outOfEquilibrium}) = flux!(f, u, edge, data, data.flux_approximation)


# The classical Scharfetter-Gummel flux scheme. This also works for space-dependent band-edge energy, but not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{scharfetter_gummel})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]

        j0                 =  params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        etak, etal         = etaFunction(u, edge, data, icc, ipsi ) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        bp, bm             = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT)
        f[icc]             = - j0 * ( bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak) )
    
    end

    return f

end

# The classical Scharfetter-Gummel flux scheme for 
# possible space-dependent DOS and band-edge energies. For these parameters the discretization scheme is modified. 
function flux!(f, u, edge, data, ::Type{scharfetter_gummel_graded})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 =  params.chargeNumbers[icc] * q * params.UT

        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]

        etak, etal         = etaFunction(u, edge, data, icc, ipsi ) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobility           = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2
        
        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT ) 
        else
            bp, bm         = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / params.UT - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc, nodek])) ) 
        end

        f[icc]             = - j0  * mobility * ( bm  * Nil * data.F[icc](etal) - bp *  Nik * data.F[icc](etak) )

    end

    return f

end

# The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{excess_chemical_potential})

    params      =   data.params
    paramsnodal =   data.paramsnodal
    
    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    

    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]
        etak, etal         = etaFunction(u, edge, data, icc, ipsi) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        Q                  = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak) )

    end

    return f

end

# The excess chemical potential flux scheme for 
# possible space-dependent DOS and band-edge energies. For these parameters the discretization scheme is modified.
function flux!(f, u, edge, data, ::Type{excess_chemical_potential_graded})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    dpsiEps     =   (params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2) * dpsi
    f[ipsi]     = - ε0 * dpsiEps
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.UT
 
        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]

        etak, etal         = etaFunction(u, edge, data, icc, ipsi ) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        bandEdgeDifference = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
        mobilityl          = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodel])
        mobilityk          = (params.mobility[icc, ireg] + paramsnodal.mobility[icc, nodek])

        if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
            Q              = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )

        else
            Q              = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc,nodek])) )

        end

        bp, bm             = fbernoulli_pm(Q)
        f[icc]             = - j0  * ( bm  * mobilityl * Nil * data.F[icc](etal) - bp * mobilityk * Nik * data.F[icc](etak) )
    
    end

    return f

end

# The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is 
# used for the regularization of the removable singularity. This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{diffusion_enhanced})

    params      =   data.params
    paramsnodal =   data.paramsnodal

    tolReg      =   1.0e-13
    
    ipsi        =   data.index_psi
    nodel       =   edge.node[2]
    nodek       =   edge.node[1]
    ireg        =   edge.region
    
    dpsi        =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]     = - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                 = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT
        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]

        etak, etal         = etaFunction(u, edge, data, icc, ipsi ) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

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
        f[icc]             = - j0 * g * (  bm * Nil * data.F[icc](etal) - bp * Nik * data.F[icc](etak))
    end

    return f

end

# The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation which arise
# when considering the generalized Scharfetter-Gummel scheme in case of Blakemore statistics.
# Hence, it should be exclusively worked with, when considering the Blakemore distribution.
# This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function flux!(f, u, edge, data, ::Type{generalized_sg})

    params        =   data.params
    paramsnodal   =   data.paramsnodal

    max_iteration =   200          # for Newton solver
    it            =   0            # number of iterations (newton)
    damp          =   0.1          # damping factor
    
    ipsi          =   data.index_psi
    nodel         =   edge.node[2]
    nodek         =   edge.node[1]
    ireg          =   edge.region
    
    dpsi          =   u[ipsi, 2] - u[ipsi, 1]
    f[ipsi]       =  - params.dielectricConstant[ireg] * ε0 * dpsi
    
    
    # k = 1 refers to left side, where as l = 2 refers to right side.
    for icc ∈ data.chargeCarrierList

        j0                  = params.chargeNumbers[icc] * q * params.mobility[icc, ireg] * params.UT

        get_DOS!(icc, edge, data)

        Nik                = data.tempDOS1[icc]
        Nil                = data.tempDOS2[icc]

        etak, etal         = etaFunction(u, edge, data, icc, ipsi ) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        bandEdgeDifference  = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

        # use Sedan flux as starting guess
        Q                   = params.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm              = fbernoulli_pm(Q)
        jInitial            = ( bm * Nik * data.F[icc](etal)  - bp * Nil * data.F[icc](etak))

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