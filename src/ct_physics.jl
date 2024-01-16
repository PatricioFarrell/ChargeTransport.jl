##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the effective DOS for interior nodes (analogously for boundary nodes and edges).
"""
function get_DOS!(icc::QType, node::VoronoiFVM.Node, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, node.region] + data.paramsnodal.densityOfStates[icc, node.index]

end

# Defining locally the effective DOS for boundary nodes.
function get_DOS!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempDOS1[icc] = data.params.bDensityOfStates[icc, bnode.region] + data.paramsnodal.densityOfStates[icc, bnode.index]

end

# Defining locally the effective DOS for edges.
function get_DOS!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempDOS1[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[1]]
    data.tempDOS2[icc] = data.params.densityOfStates[icc, edge.region] + data.paramsnodal.densityOfStates[icc, edge.node[2]]

    return data.tempDOS1, data.tempDOS2
end

# Calculate the DOS on a given interior region.
function get_DOS(icc::QType, ireg::Int, ctsys)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    return data.params.densityOfStates[icc, ireg] .+ view(data.paramsnodal.densityOfStates[icc, :], subgrid(grid, [ireg])) # view nodal dependent DOS on respective grid
end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Defining locally the band-edge energy for interior nodes (analougesly for boundary nodes and edges).
"""
function get_BEE!(icc::QType, node::VoronoiFVM.Node, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, node.region] + data.paramsnodal.bandEdgeEnergy[icc, node.index]

end

# Defining locally the band-edge energy for boundary nodes.
function get_BEE!(icc::QType, bnode::VoronoiFVM.BNode, data)

    data.tempBEE1[icc] = data.params.bBandEdgeEnergy[icc, bnode.region] + data.paramsnodal.bandEdgeEnergy[icc, bnode.index]

end

# Defining locally the band-edge energy for edges.
function get_BEE!(icc::QType, edge::VoronoiFVM.Edge, data)

    data.tempBEE1[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[1]]
    data.tempBEE2[icc] = data.params.bandEdgeEnergy[icc, edge.region] + data.paramsnodal.bandEdgeEnergy[icc, edge.node[2]]

    return data.tempBEE1, data.tempBEE2

end

# Calculate the band-edge energy on a given interior region.
function get_BEE(icc::QType, ireg::Int, ctsys)

    grid = ctsys.fvmsys.grid
    data = ctsys.fvmsys.physics.data

    return data.params.bandEdgeEnergy[icc, ireg] .+ view(data.paramsnodal.bandEdgeEnergy[icc, :], subgrid(grid, [ireg]))

end
##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

The argument of the statistics function for interior nodes.
"""
function etaFunction!(u, node::VoronoiFVM.Node, data, icc)

    get_BEE!(icc, node::VoronoiFVM.Node, data)

    E  = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[data.index_psi]) + E / q )

end


### DA: DELETE THISSSSS
"""
$(TYPEDSIGNATURES)

The argument of the distribution function for floats.
"""
function etaFunction(u, data, node, region, icc, in_region::Bool)

    if in_region == true
        E  = data.params.bandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    elseif in_region == false
        E  = data.params.bBandEdgeEnergy[icc, region] + data.paramsnodal.bandEdgeEnergy[icc, node]
    end

    return etaFunction(u[data.index_psi], u[icc], data.params.UT, E, data.params.chargeNumbers[icc])
end


"""
$(TYPEDSIGNATURES)

The argument of the statistics function for boundary nodes.
"""
function etaFunction!(u, bnode::VoronoiFVM.BNode, data, icc) # bnode.index refers to index in overall mesh

    get_BEE!(icc, bnode::VoronoiFVM.BNode, data)
    E  = data.tempBEE1[icc]

    return data.params.chargeNumbers[icc] / data.params.UT * ( (u[icc] - u[data.index_psi]) + E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the statistics function for edges.
"""

function etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    get_BEE!(icc, edge::VoronoiFVM.Edge, data)

    E1 = data.tempBEE1[icc];  E2 = data.tempBEE2[icc]

    return etaFunction(u[data.index_psi, 1], u[icc, 1], data.params.UT, E1, data.params.chargeNumbers[icc]), etaFunction(u[data.index_psi, 2], u[icc, 2], data.params.UT, E2, data.params.chargeNumbers[icc])
end

"""
$(TYPEDSIGNATURES)

The argument of the statistics function for a given solution on a given interior region.
"""
function etaFunction(sol, ireg::Int, ctsys, icc::QType)

    grid   = ctsys.fvmsys.grid
    data   = ctsys.fvmsys.physics.data

    Ecc    = get_BEE(icc, ireg, ctsys)
    # view solution on respective grid
    solcc  = view(sol[icc, :], subgrid(grid, [ireg]))
    solpsi = view(sol[data.index_psi, :], subgrid(grid, [ireg]))

    return data.params.chargeNumbers[icc] ./ data.params.UT .* ( (solcc .- solpsi) .+ Ecc ./ q )
end


"""
$(TYPEDSIGNATURES)

The argument of the statistics function for given ``\\varphi_\\alpha``
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

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for interior nodes.

"""
function get_density!(u, node::VoronoiFVM.Node, data, icc)

    get_DOS!(icc, node, data)

    Ncc   = data.tempDOS1[icc]
    eta   = etaFunction!(u, node, data, icc) # calls etaFunction!(u,node::VoronoiFVM.Node,data,icc)

    return Ncc * data.F[icc](eta)

end

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for interior nodes.

"""
function get_density!(u, bnode::VoronoiFVM.BNode, data, icc)

    get_DOS!(icc, bnode, data)

    Ncc   = data.tempDOS1[icc]
    eta   = etaFunction!(u, bnode, data, icc) # calls etaFunction!(u,node::VoronoiFVM.BNode,data,icc)

    return Ncc * data.F[icc](eta)

end


"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for edges.

"""
function get_density!(u, edge::VoronoiFVM.Edge, data, icc)

    get_DOS!(icc, edge, data)

    Ncc1       = data.tempDOS1[icc]
    Ncc2       = data.tempDOS2[icc]

    eta1, eta2 = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    return Ncc1 * data.F[icc](eta1), Ncc2 * data.F[icc](eta2)

end


"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities for given interior region corresponding
to a homogeneous set of parameters.

"""
function get_density(sol, ireg::Int, ctsys, icc::QType)

    data = ctsys.fvmsys.physics.data

    Ncc  = get_DOS(icc, ireg, ctsys)
    eta  = etaFunction(sol, ireg, ctsys, icc)

    return Ncc .* data.F[icc].(eta)

end
###########################################################
###########################################################

function emptyFunction()
end

"""
Function in case of an applied voltage equal to zero at one boundary.
"""
zeroVoltage(t) = 0.0

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master breaction! function. This is the function which enters VoronoiFVM and hands over
for each boundary the chosen boundary model.

"""
breaction!(f, u, bnode, data) = breaction!(f, u, bnode, data, data.boundaryType[bnode.region])

#################################################################################################

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
function breaction!(f, u, bnode, data, ::Type{OhmicContact})

    params      = data.params
    paramsnodal = data.paramsnodal

    ipsi        = data.index_psi

    # electrons and holes entering right hand-side for BC of ipsi
    for icc ∈ data.electricCarrierList         # Array{Int64, 1}

        icc     = data.chargeCarrierList[icc]  # Array{QType, 1}
        ncc     = get_density!(u, bnode, data, icc)

        # subtract doping
        f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )
        # add charge carrier
        f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

    end

    # if ionic carriers are present
    for iicc ∈ data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if bnode.cellregions[1] ∈ iicc.regions    # bnode.cellregions = [bnode.region, 0] for outer boundary.
            icc     = iicc.ionicCarrier           # species number chosen by user
            icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            ncc     = get_density!(u, bnode, data, icc)

            # subtract doping
            f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )
            # add charge carrier
            f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

        end

    end

    # # if trap carriers are present
    # for iicc ∈ data.trapCarrierList # ∈ Array{TrapCarrier, 1}
    #     # add trap carriers only in defined regions (otherwise get NaN error)
    #     if bnode.cellregions[1] ∈ iicc.regions
    #         icc     = iicc.trapCarrier           # species number chosen by user
    #         icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

    #         ncc     = get_density!(u, bnode, data, icc)

    #         # subtract doping
    #         f[ipsi] = f[ipsi] - params.chargeNumbers[icc] * ( params.bDoping[icc, bnode.region] )
    #         # add charge carrier
    #         f[ipsi] = f[ipsi] + params.chargeNumbers[icc] * ncc

    #     end

    # end

    f[ipsi] = f[ipsi] - paramsnodal.doping[bnode.index]
    f[ipsi] = - data.λ1 * 1 / tiny_penalty_value *  q * f[ipsi]

    # electrons and holes boundary condition
    iphin    = data.bulkRecombination.iphin # integer index of φ_n
    iphip    = data.bulkRecombination.iphip # integer index of φ_p

    Δu       = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    boundary_dirichlet!(f, u, bnode, species = iphin, region=bnode.region, value=Δu)
    boundary_dirichlet!(f, u, bnode, species = iphip, region=bnode.region, value=Δu)

end

"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions. For the electrostatic potential we assume

``\\psi = - \\phi_S/q + U, ``

where  ``\\phi_S`` corresponds to a given value (non-negative Schottky barrier) and ``U`` to the applied voltage. The quantitity ``\\phi_S`` needs to be specified in the main file.
For eletrons and holes we assume the following

``f[n_\\alpha]  =  z_\\alpha q v_\\alpha (n_\\alpha - n_{\\alpha, 0})``,

where ``v_{\\alpha}`` can be treated as a surface recombination mechanism and is given. The parameter
``n_{\\alpha, 0}`` is the equilibrium density of the charge carrier ``\\alpha`` and can be
calculated via

``n_{\\alpha, 0}= N_\\alpha \\mathcal{F}_\\alpha \\Bigl( - z_\\alpha/ U_T (E_c - E_\\alpha) - \\phi_S) / q \\Bigr). ``

"""

function breaction!(f, u, bnode, data, ::Type{SchottkyContact})

    params      = data.params
    paramsnodal = data.paramsnodal
    ipsi        = data.index_psi
    iphin       = data.bulkRecombination.iphin
    Ec          = params.bBandEdgeEnergy[iphin, bnode.region]

    for icc ∈ data.electricCarrierList       # Array{Int64, 1}

        icc    = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used

        get_DOS!(icc, bnode, data);  get_BEE!(icc, bnode, data)
        Ni     =   data.tempDOS1[icc]
        Ei     =   data.tempBEE1[icc]
        etaFix = - params.chargeNumbers[icc] / params.UT * (  ( (Ec - Ei) - params.SchottkyBarrier[bnode.region] ) / q  )

        ncc    = get_density!(u, bnode, data, icc)

        f[icc] =   params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * ( ncc - Ni * data.F[icc](etaFix) )

    end

    # function evaluation causes allocation!!!
    Δu         = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    ipsiIndex  = length(data.chargeCarrierList) + 1 # This is necessary, since passing something other than an Integer in boundary_dirichlet!() causes allocations
    boundary_dirichlet!(f, u, bnode, species=ipsiIndex, region=bnode.region, value=(- (params.SchottkyBarrier[bnode.region]  - Ec)/q) + Δu)

end


"""
$(TYPEDSIGNATURES)
Creates Schottky boundary conditions with additional lowering which are modelled as

`` \\psi = - \\phi_S/q  + \\sqrt{ -\\frac{ q  \\nabla_{\\boldsymbol{\\nu}} \\psi_\\mathrm{R}}{4\\pi \\varepsilon_\\mathrm{i}}} + U``,

where `` \\psi_\\mathrm{R}`` denotes the electric potential with standard Schottky contacts and the same space charge density as `` \\psi`` and where ``\\varepsilon_\\mathrm{i}}}`` corresponds to the image force dielectric constant.

To solve for this additional boundary conditions the projected gradient ``\\nabla_{\\boldsymbol{\\nu}} \\psi_\\mathrm{R} `` is stored within a boundary species and calculated in the method generic_operator!().
"""
function breaction!(f, u, bnode, data, ::Type{SchottkyBarrierLowering})

    params       = data.params
    ipsi         = data.index_psi
    iphin        = data.bulkRecombination.iphin
    Ec           = params.bBandEdgeEnergy[iphin, bnode.region]
    ipsiStandard = data.barrierLoweringInfo.ipsiStandard
    ipsiGrad     = data.barrierLoweringInfo.ipsiGrad

    if data.calculationType == OutOfEquilibrium
        for icc ∈ data.electricCarrierList       # Array{Int64, 1}

            icc    = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used
            ncc    = get_density!(u, bnode, data, icc)

            f[icc] = params.chargeNumbers[icc] * q *  params.bVelocity[icc, bnode.region] * ( ncc - params.bDensityEQ[icc, bnode.region] )

        end
    end

    # function evaluation causes allocation!!!
    Δu         = params.contactVoltage[bnode.region] + data.contactVoltageFunction[bnode.region](bnode.time)

    if u[ipsiGrad] < 0
        PsiS = sqrt( - q/(4 * pi * params.dielectricConstantImageForce[bnode.cellregions[1]]) * u[ipsiGrad]) #bnode.cellregions[1] ∈ iicc.regions    # bnode.cellregions = [bnode.region, 0] for outer boundary.
    else
        PsiS = 0.0
    end

    f[ipsiStandard] = 1/tiny_penalty_value * (u[ipsi]  + (params.SchottkyBarrier[bnode.region] - Ec)/q - PsiS - Δu)
    f[ipsi]         = 1/tiny_penalty_value * (u[ipsiStandard] + (params.SchottkyBarrier[bnode.region] - Ec)/q        - Δu)

end


# This breaction! function is chosen when no interface model is chosen.
breaction!(f, u, bnode, data, ::Type{InterfaceNone}) = emptyFunction()


function breaction!(f, u, bnode, data, ::Type{InterfaceRecombination})

    params = data.params

    if data.calculationType == InEquilibrium
        return
    end

    # indices (∈ IN) of electron and hole quasi Fermi potentials specified by user (passed through recombination)
    iphin           = data.bulkRecombination.iphin # integer index of φ_n
    iphip           = data.bulkRecombination.iphip # integer index of φ_p

    n               = get_density!(u, bnode, data, iphin)
    p               = get_density!(u, bnode, data, iphip)

    exponentialTerm = exp((q * u[iphin] - q  * u[iphip] ) / (kB * params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    kernelSRH = 1.0 / ( 1.0/params.recombinationSRHvelocity[iphip, bnode.region] * (n + params.bRecombinationSRHTrapDensity[iphin, bnode.region]) + 1.0/params.recombinationSRHvelocity[iphin, bnode.region] * (p + params.bRecombinationSRHTrapDensity[iphip, bnode.region] ) )

    for icc ∈ data.electricCarrierList
        icc = data.chargeCarrierList[icc]
        f[icc] =  q * params.chargeNumbers[icc] * kernelSRH *  excessDensTerm
    end


end

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Generic operator to save the projected gradient of electric potential
(for system with standard Schotty contacts). Note that this currently
only working in one dimension!

"""
function generic_operator!(f, u, fvmsys)

    f                  .=0

    coord               = fvmsys.grid[Coordinates]
    n                   = length(coord)
    barrierLoweringInfo = fvmsys.physics.data.barrierLoweringInfo
    ipsiStandard        = barrierLoweringInfo.ipsiStandard
    ipsiGrad            = barrierLoweringInfo.ipsiGrad

    idx                 = barrierLoweringInfo.idx

    f[idx[ipsiGrad, 1]] = u[idx[ipsiGrad, 1]] - (-1) * ( u[idx[ipsiStandard, 2]] - u[idx[ipsiStandard, 1]]   ) / (coord[2] - coord[1]  )
    f[idx[ipsiGrad, n]] = u[idx[ipsiGrad, n]] -        ( u[idx[ipsiStandard, n]] - u[idx[ipsiStandard, n-1]] ) / (coord[n] - coord[n-1])


end


##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bstorage! function. This is the function which enters VoronoiFVM and hands over
for each boundary the time-dependent part of the chosen boundary model.

"""
bstorage!(f, u, bnode, data) = bstorage!(f, u, bnode, data, data.modelType)

bstorage!(f, u, bnode, data, ::Type{Stationary})  = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{Transient}) = bstorage!(f, u, bnode, data, data.boundaryType[bnode.region])


bstorage!(f, u, bnode, data, ::Type{InterfaceNone}) = emptyFunction()

bstorage!(f, u, bnode, data, ::Type{InterfaceRecombination}) = emptyFunction()

# No bstorage! is used, if an ohmic and schottky contact model is chosen.
bstorage!(f, u, bnode, data, ::OuterBoundaryModelType) =  emptyFunction()

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master bflux! function. This is the function which enters VoronoiFVM and hands over
for each boundary the flux within the boundary.

"""
bflux!(f, u, bedge, data) = emptyFunction()



##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master reaction! function. This is the function which enters VoronoiFVM and hands over
reaction terms for concrete calculation type and bulk recombination model.

"""
reaction!(f, u, node, data) = reaction!(f, u, node, data, data.calculationType)

"""
$(TYPEDSIGNATURES)
Reaction in case of equilibrium, i.e. no generation and recombination is considered.
"""
function reaction!(f, u, node, data, ::Type{InEquilibrium})

    ipsi = data.index_psi
    # RHS of Poisson
    RHSPoisson!(f, u, node, data, ipsi)
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        RHSPoisson!(f, u, node, data, ipsiStandard)
    end

    # zero reaction term for all icc (stability purpose)
    for icc ∈ data.electricCarrierList # Array{Int64, 1}
        icc = data.chargeCarrierList[icc] # Array{QType 1}
        f[icc] = u[icc]
    end

    for iicc ∈ data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions
            icc     = iicc.ionicCarrier           # species number chosen by user
            icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            f[icc] = u[icc]
        end
    end

    for iicc ∈ data.trapCarrierList # ∈ Array{TrapCarrier, 1}
        # add trap carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions
            icc     = iicc.trapCarrier            # species number chosen by user
            icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

            f[icc] = u[icc]
        end
    end

end


function addRecombination!(f, u, node, data, ::SRHWithoutTrapsType)

    params = data.params
    ireg   = node.region

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin  = data.bulkRecombination.iphin
    iphip  = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin  = data.chargeCarrierList[iphin]
    iphip  = data.chargeCarrierList[iphip]

    n      = get_density!(u, node, data, iphin)
    p      = get_density!(u, node, data, iphip)

    taun   = params.recombinationSRHLifetime[iphin, ireg]
    n0     = params.recombinationSRHTrapDensity[iphin, ireg]
    taup   = params.recombinationSRHLifetime[iphip, ireg]
    p0     = params.recombinationSRHTrapDensity[iphip, ireg]

    exponentialTerm = exp((q * u[iphin] - q * u[iphip]) / (kB * data.params.temperature))
    excessDensTerm  = n * p * (1.0 - exponentialTerm)

    # calculate recombination kernel. If user adjusted Auger, radiative or SRH recombination,
    # they are set to 0. Hence, adding them here, has no influence since we simply add by 0.0.
    kernelRad   = params.recombinationRadiative[ireg]
    kernelAuger = (params.recombinationAuger[iphin, ireg] * n + params.recombinationAuger[iphip, ireg] * p)
    kernelSRH   = params.prefactor_SRH / ( taup * (n + n0) + taun * (p + p0) )
    kernel      = kernelRad + kernelAuger + kernelSRH
    ###########################################################
    ####       right-hand side of continuity equations     ####
    ####       for φ_n and φ_p (bipolar reaction)          ####
    ###########################################################
    f[iphin] = q * params.chargeNumbers[iphin] * kernel * excessDensTerm
    f[iphip] = q * params.chargeNumbers[iphip] * kernel * excessDensTerm

end


function addRecombination!(f, u, node, data, ::SRHWithTrapsType)

    params = data.params
    ireg   = node.region

    # indices (∈IN) used by user
    iphin  = data.bulkRecombination.iphin
    iphip  = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin  = data.chargeCarrierList[iphin]
    iphip  = data.chargeCarrierList[iphip]

    n      = get_density!(u, node, data, iphin)
    p      = get_density!(u, node, data, iphip)

    taun   = params.recombinationSRHLifetime[iphin, ireg]
    n0     = params.recombinationSRHTrapDensity[iphin, ireg]
    taup   = params.recombinationSRHLifetime[iphip, ireg]
    p0     = params.recombinationSRHTrapDensity[iphip, ireg]

    for iicc ∈ data.trapCarrierList
        # add trap carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions

            icc   = iicc.trapCarrier
            itrap = data.chargeCarrierList[icc]

            Nt    = params.densityOfStates[itrap, ireg]
            t     = get_density!(u, node, data, itrap)

            # Rn, Rp agree up to sign with *On the Shockley-Read-Hall Model: Generation-Recombination
            # in Semiconductors* in SIAM Journal on Applied Mathematics, Vol. 67, No. 4 (2007), pp. 1183-1201.
            # The sign is chosen according to *Supporting Information: Consistent Device Simulation Model
            # Describing Perovskite Solar Cells in Steady-State, Transient and Frequency Domain* in ACS (2018)
            if params.chargeNumbers[itrap] == -1
                Rn = 1 / taun * (n * (1-t/Nt) - n0 * t/Nt)
                Rp = 1 / taup * (p * t/Nt     - p0 * (1-t/Nt))
            elseif params.chargeNumbers[itrap] == 1
                Rn = 1 / taun * (n * t/Nt     - n0 * (1-t/Nt))
                Rp = 1 / taup * (p * (1-t/Nt) - p0 * t/Nt)
            end

            f[itrap] = q * params.chargeNumbers[itrap] * (Rp-Rn)
            f[iphin] = q * params.chargeNumbers[iphin] * Rn
            f[iphip] = q * params.chargeNumbers[iphip] * Rp
        end

    end

end

function addGeneration!(f, u, node, data)

    generationTerm = generation(data, node.region, node.coord[node.index], data.generationModel)

    for icc ∈ data.electricCarrierList
        icc    = data.chargeCarrierList[icc] # based on user index and regularity of solution quantities or integers are used and depicted here
        f[icc] = f[icc] - q * data.params.chargeNumbers[icc] * generationTerm
    end


end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of Poisson equation, i.e. which builds
the space charge density.
"""
function RHSPoisson!(f, u, node, data, ipsi)

    ###########################################################
    ####         right-hand side of nonlinear Poisson      ####
    ####         equation (space charge density)           ####
    ###########################################################

    # electrons and holes entering right hand-side of Poisson in each layer
    for icc ∈ data.electricCarrierList          # Array{Int64, 1}

        icc     = data.chargeCarrierList[icc]   # Array{QType, 1}
        ncc     = get_density!(u, node, data, icc)

        f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * ( data.params.doping[icc, node.region] )  # subtract doping
        f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier

    end

    for iicc ∈ data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # add ionic carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions

            icc     = iicc.ionicCarrier           # species number chosen by user
            icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})
            ncc     = get_density!(u, node, data, icc)

            f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * ( data.params.doping[icc, node.region] )  # subtract doping
            f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier
        end
    end

    for iicc ∈ data.trapCarrierList # ∈ Array{TrapCarrier, 1}
        # add trap carriers only in defined regions (otherwise get NaN error)
        if node.region ∈ iicc.regions

            icc     = iicc.trapCarrier            # species number chosen by user
            icc     = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})
            ncc     = get_density!(u, node, data, icc)

            f[ipsi] = f[ipsi] - data.params.chargeNumbers[icc] * ( data.params.doping[icc, node.region] )  # subtract doping
            f[ipsi] = f[ipsi] + data.params.chargeNumbers[icc] * ncc   # add charge carrier
        end
    end

    f[ipsi] = f[ipsi] - data.paramsnodal.doping[node.index]

    f[ipsi] = - q * data.λ1 * f[ipsi]

    ## This is the trap density for the stationary case without traps as own charge carrier
    addTrapDensity!(f, u, node, data)

end

"""
$(TYPEDSIGNATURES)
Function which builds right-hand side of electric charge carriers.
"""
function RHSContinuityEquations!(f, u, node, data)

    # dependent on user information concerncing recombination
    addRecombination!(f, u, node, data, data.bulkRecombination.bulk_recomb_SRH)
    # dependent on user information concerncing generation
    addGeneration!(f, u, node, data)

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
function reaction!(f, u, node, data, ::Type{OutOfEquilibrium})

    ipsi = data.index_psi
    RHSPoisson!(f, u, node, data, ipsi)               # RHS of Poisson
    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn # additional Poisson in case of Barrier Lowering
        ipsiStandard = data.barrierLoweringInfo.ipsiStandard
        RHSPoisson!(f, u, node, data, ipsiStandard)
    end

    # First, set RHS to zero for all icc
    for icc ∈ eachindex(data.chargeCarrierList) # Array{Int61, 1}
        icc     = data.chargeCarrierList[icc]   # Array{QType, 1}
        f[icc]  = 0.0
    end

    # Then, add RHS of continuity equations based on user information
    RHSContinuityEquations!(f, u, node, data) # RHS of Charge Carriers with special treatment of recombination

end

# Function which adds additional trap density to right-hand side of Poisson equation
# without modeling traps as own charge carrier.
# Note that this one may be deleted in future version.
addTrapDensity!(f, u, node, data) = addTrapDensity!(f, u, node, data, data.bulkRecombination.SRH_2species_trap)


function addTrapDensity!(f, u, node, data, ::SRHWithoutTrapsType)
    return
end

function addTrapDensity!(f, u, node, data, ::Type{SRH2SpeciesPresentTrapDens})

    params = data.params
    ireg   = node.region
    ipsi   = data.index_psi

    # indices (∈ IN) of electron and hole quasi Fermi potentials used by user (passed through recombination)
    iphin  = data.bulkRecombination.iphin
    iphip  = data.bulkRecombination.iphip

    # based on user index and regularity of solution quantities or integers are used and depicted here
    iphin  = data.chargeCarrierList[iphin]
    iphip  = data.chargeCarrierList[iphip]

    n      = get_density!(u, node, data, iphin)
    p      = get_density!(u, node, data, iphip)

    n0     = params.recombinationSRHTrapDensity[iphin, ireg]
    p0     = params.recombinationSRHTrapDensity[iphip, ireg]
    taun   = params.recombinationSRHLifetime[iphin, ireg]
    taup   = params.recombinationSRHLifetime[iphip, ireg]
    zt     = data.AuxTrapValues.zt
    Nt     = data.AuxTrapValues.Nt[ireg]


    # add trap density
    if zt==1
        f[ipsi] = f[ipsi] - q * zt * data.λ1 * Nt * ( 1 - (taun*p0 + taup*n) / (taun*(p0+p) + taup*(n0+n)) )
    elseif zt==-1
        f[ipsi] = f[ipsi] + q * zt * data.λ1 * Nt * ( (taun*p0 + taup*n) / (taun*(p0+p) + taup*(n0+n)) )
    end

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
function generation(data, ireg, node, ::Type{GenerationUniform})

    return data.λ2 * data.params.generationUniform[ireg]
end


# The generation rate ``G``, which occurs in the right-hand side of the
# continuity equations obeying the Beer-Lambert law.
# only works in 1D till now; adjust node, when multidimensions
function generation(data, ireg, node, ::Type{GenerationBeerLambert})

    params = data.params

    return data.λ2 * params.generationIncidentPhotonFlux[ireg] * params.generationAbsorption[ireg] * exp( - params.invertedIllumination * params.generationAbsorption[ireg] * (node - params.generationPeak))

end

generation(data, ireg, node, ::Type{GenerationNone}) = 0.0

##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)
Master storage! function. This is the function which enters VoronoiFVM and hands over
a storage term, if we consider transient problem.
"""
storage!(f, u, node, data) = storage!(f, u, node, data, data.modelType)

storage!(f, u, node, data, ::Type{Stationary})  = emptyFunction()


storage!(f, u, node, data, ::Type{Transient}) = storage!(f, u, node, data, data.calculationType)

storage!(f, u, node, data, ::Type{InEquilibrium}) = emptyFunction()
"""
$(TYPEDSIGNATURES)
The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have
``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha``
and for the electrostatic potential
``f[ψ] = 0``.
"""
function storage!(f, u, node, data, ::Type{OutOfEquilibrium})

    params = data.params
    ipsi   = data.index_psi

    for icc ∈ data.electricCarrierList       # Array{Int64, 1}

        icc    = data.chargeCarrierList[icc] # get correct index in chargeCarrierList
        ncc    = get_density!(u, node, data, icc)
        f[icc] = q * params.chargeNumbers[icc] * ncc

    end

    for iicc ∈ data.ionicCarrierList # ∈ Array{IonicCarrier, 1}
        # Here we do not need to check, if carrier is present in a specific region.
        # This is directly handled by VoronoiFVM.
        icc    = iicc.ionicCarrier           # species number chosen by user
        icc    = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

        ncc    = get_density!(u, node, data, icc)
        f[icc] = q * params.chargeNumbers[icc] * ncc
    end

    for iicc ∈ data.trapCarrierList # ∈ Array{TrapCarrier, 1}
        # Here we do not need to check, if carrier is present in a specific region.
        # This is directly handled by VoronoiFVM.
        icc    = iicc.trapCarrier            # species number chosen by user
        icc    = data.chargeCarrierList[icc] # find correct index within chargeCarrierList (Array{QType, 1})

        ncc    = get_density!(u, node, data, icc)

        f[icc] = q * params.chargeNumbers[icc] * ncc
    end


    f[ipsi] = 0.0

end

##########################################################
##########################################################
"""
$(TYPEDSIGNATURES)
Master flux functions which enters VoronoiFVM. Flux discretization scheme is chosen in two steps. First, we need
to see, if we are in or out of equilibrium. If, InEquilibrium, then
no flux is passed. If outOfEquilibrium, we choose the flux approximation
which the user chose for each charge carrier. For the displacement flux we use a finite difference approach.
"""
flux!(f, u, edge, data) = flux!(f, u, edge, data, data.calculationType)


# Finite difference discretization of the displacement flux.
function displacementFlux!(f, u, edge, data)

    params      =   data.params
    paramsnodal =   data.paramsnodal

    ipsi        =   data.index_psi
    nodel       =   edge.node[2]   # left node
    nodek       =   edge.node[1]   # right node
    ireg        =   edge.region

    dpsi        =   u[ipsi, 2] - u[ipsi, 1]

    dielConst   =   params.dielectricConstant[ireg]  + (paramsnodal.dielectricConstant[nodel] + paramsnodal.dielectricConstant[nodek])/2
    f[ipsi]     = - dielConst * dpsi

    if data.barrierLoweringInfo.BarrierLoweringOn == BarrierLoweringOn # additional Poisson in case of Barrier Lowering
        ipsiStandard    = data.barrierLoweringInfo.ipsiStandard
        f[ipsiStandard] = - dielConst * (u[ipsiStandard, 2] - u[ipsiStandard, 1])
    end

end


function flux!(f, u, edge, data, ::Type{InEquilibrium})
    ## discretization of the displacement flux (LHS of Poisson equation)
    displacementFlux!(f, u, edge, data)
end

function flux!(f, u, edge, data, ::Type{OutOfEquilibrium})

    ## discretization of the displacement flux (LHS of Poisson equation)
    displacementFlux!(f, u, edge, data)

    for icc ∈ data.electricCarrierList   # correct index of electric carriers of Type Int64
        chargeCarrierFlux!(f, u, edge, data, icc, data.fluxApproximation[icc])
    end

    for icc ∈ data.ionicCarrierList
        icc     = icc.ionicCarrier       # correct index number chosen by user of Type Int64
        chargeCarrierFlux!(f, u, edge, data, icc, data.fluxApproximation[icc])
    end

    for icc ∈ data.trapCarrierList
        icc     = icc.trapCarrier       # correct index number chosen by user of Type Int64
        chargeCarrierFlux!(f, u, edge, data, icc, data.fluxApproximation[icc])
    end

end


# The classical Scharfetter-Gummel flux scheme. This also works for space-dependent
# band-edge energy, but not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ScharfetterGummel})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region
    j0           = params.UT * params.mobility[icc, ireg]

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    bp, bm       = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT)
    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp * ncck )

end

# The classical Scharfetter-Gummel flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the
# discretization scheme is modified.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ScharfetterGummelGraded})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region

    mobility     = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2
    j0           = params.UT * mobility

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
        bp, bm   = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT )
    else
        bp, bm   = fbernoulli_pm( params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / params.UT - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc, nodek])) )
    end

    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp *  ncck )

end

# The excess chemical potential flux discretization scheme. This also works for space-dependent band-edge energy, but
# not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ExcessChemicalPotential})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region
    j0           = params.UT * params.mobility[icc, ireg]

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal   = etaFunction!(u, edge, data, icc)

    Q            = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff /q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
    bp, bm       = fbernoulli_pm(Q)

    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp * ncck )

end

# The excess chemical potential flux scheme for
# possible space-dependent DOS and band-edge energies. For these parameters the discretization scheme is modified.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{ExcessChemicalPotentialGraded})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region

    mobility     = params.mobility[icc, ireg] + (paramsnodal.mobility[icc, nodel] + paramsnodal.mobility[icc, nodek])/2
    j0           = params.UT * mobility

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal   = etaFunction!(u, edge, data, icc)

    if paramsnodal.densityOfStates[icc, nodel] ≈ 0.0 || paramsnodal.densityOfStates[icc, nodek] ≈ 0.0
        Q        = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
    else
        Q        = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(paramsnodal.densityOfStates[icc, nodel]) -log(paramsnodal.densityOfStates[icc,nodek])) )
    end

    bp, bm       = fbernoulli_pm(Q)
    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp * ncck )

end

# The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme
# is used for the regularization of the removable singularity. This also works for
# space-dependent band-edge energy, but not for space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{DiffusionEnhanced})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal   = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    if ( log(data.F[icc](etal)) - log(data.F[icc](etak)) ) ≈ 0.0 # regularization idea coming from Pietra-Jüngel scheme
        gk = exp(etak) / data.F[icc](etak)
        gl = exp(etal) / data.F[icc](etal)
        g  = 0.5 * ( gk + gl )
    else
        g  = (etal - etak ) / ( log(data.F[icc](etal)) - log(data.F[icc](etak)) )
    end

    j0           = params.UT * params.mobility[icc, ireg] * g

    bp, bm       = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / (params.UT * g))
    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp * ncck)

end

# The diffusion enhanced scheme by Bessemoulin-Chatard for fluxes not based on a nonlinear diffusion
# but on a modified drift.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{DiffusionEnhancedModifiedDrift})

    params       = data.params
    paramsnodal  = data.paramsnodal

    icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]

    etak, etal   = etaFunction!(u, edge, data, icc) # calls etaFunction!(u, edge::VoronoiFVM.Edge, data, icc)

    if ( log(data.F[icc](etal)) - log(data.F[icc](etak)) ) ≈ 0.0 # regularization idea coming from Pietra-Jüngel scheme
        gk = exp(etak) / data.F[icc](etak)
        gl = exp(etal) / data.F[icc](etal)
        g  = 0.5 * ( gk + gl )
    else
        g  = (etal - etak ) / ( log(data.F[icc](etal)) - log(data.F[icc](etak)) )
    end

    j0           = params.UT * params.mobility[icc, ireg]

    bp, bm       = fbernoulli_pm(params.chargeNumbers[icc] * (dpsi - bandEdgeDiff / q) / (params.UT * g))
    ncck, nccl   = get_density!(u, edge, data, icc)

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * ( bm * nccl - bp * ncck)

end


# # The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation
# # which arise when considering the generalized Scharfetter-Gummel scheme in case of Blakemore
# # statistics. Hence, it should be exclusively worked with, when considering the Blakemore
# # statistics. This also works for space-dependent band-edge energy, but not for
# # space-dependent effective DOS.
function chargeCarrierFlux!(f, u, edge, data, icc, ::Type{GeneralizedSG})

    max_iter     = 300          # for Newton solver
    it           = 0            # number of iterations (newton)
    damp         = 0.1          # damping factor

    params       = data.params
    paramsnodal  = data.paramsnodal

    # DA: we get issues with allocations, when allowing non Integer icc.
    #icc          = data.chargeCarrierList[icc]
    ipsi         = data.index_psi
    nodek        = edge.node[1]   # left node
    nodel        = edge.node[2]   # right node
    ireg         = edge.region
    j0           = params.UT * params.mobility[icc, ireg]

    dpsi         = u[ipsi, 2] - u[ipsi, 1]
    bandEdgeDiff = paramsnodal.bandEdgeEnergy[icc, nodel] - paramsnodal.bandEdgeEnergy[icc, nodek]
    etak, etal   = etaFunction!(u, edge, data, icc)

    # use Sedan flux as starting guess
    Q            = params.chargeNumbers[icc]*( (dpsi - bandEdgeDiff/q) /params.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
    bp, bm       = fbernoulli_pm(Q)
    ncck, nccl   = get_density!(u, edge, data, icc)

    jInitial     = ( bm * nccl - bp * ncck)

    implicitEq(j::Real) = (fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDiff/q)) /params.UT + params.γ * j )[2] * exp(etal) - fbernoulli_pm(params.chargeNumbers[icc] * ((dpsi - bandEdgeDiff/q) /params.UT) - params.γ * j)[1] * exp(etak)) - j

    delta        = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
    oldup        = 1.0
    while (it < max_iter)
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
        oldup = value(update)

        it       = it + 1
        damp     = min(damp * 1.2, 1.0)
    end

    f[icc]       = - params.chargeNumbers[icc] * q * j0 * jInitial

end

# ##########################################################
# ##########################################################