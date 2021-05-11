"""
$(TYPEDEF)

A struct holding the physical data for a drift-diffusion simulation of a semiconductor device.
If there are ``N`` number of species ``\\alpha`` with ``\\alpha = 1, ..., N`` it is assumed that the first ``N-1`` species
correspond to the charge carriers with a density ``n_\\alpha``. The final species corresponds to the electrostatic potential ``\\psi``.

$(TYPEDFIELDS)
"""
mutable struct ChargeTransportData <: VoronoiFVM.AbstractData

    # integer numbers

    """
    number of nodes used for the disretization of the domain ``\\mathbf{\\Omega}``.
    """
    numberOfNodes               ::  Int64

    """
    number of regions ``\\mathbf{\\Omega}_k`` within the domain ``\\mathbf{\\Omega}``.
    """
    numberOfRegions             ::  Int64

    """
    number of boundary regions ``(\\partial \\mathbf{\\Omega})_k`` such that 
    `` \\partial \\mathbf{\\Omega} = \\cup_k (\\partial \\mathbf{\\Omega})_k``.
    """
    numberOfBoundaryRegions     ::  Int64

    """
    number of carriers, which corresponds to the number of species ``-1``.
    """
    numberOfCarriers            ::  Int64

    numberOfInterfaceSpecies    ::  Int64


    # real numbers

    """
    A given constant temperature.
    """
    temperature                 ::  Float64

    """
    The thermal voltage, which reads  ``U_T = k_B T / q``.
    """
    UT                          ::  Float64

    """
    A reference energy, which is only used for numerical computations.
    """
    # DA: I think that, when using Schottky contacts we need to set the reference
    # energy to the Fermi Level to get a physical reasonable electrostatic potential,
    # but I am not sure yet and did not test it well ...
    Eref                        ::  Float64

    """
    The parameter of the Blakemore statistics.
    """
    γ                           ::  Float64

    """
    An embedding parameter used to solve the nonlinear Poisson problem, which results
    in the case of thermodynamic equilibrium and electrocharge neutrality.
    """
    λ1                          ::  Float64

    """
    An embedding parameter for turning the generation rate ``G`` on.
    """
    λ2                          ::  Float64

    """
    An embedding parameter for electrochemical reaction.
    """
    λ3                          ::  Float64

    """
    Prefactor of electro-chemical reaction of internal boundary conditions.
    """
    r0                          ::  Float64


    # booleans

    """
    A boolean, which is true in case of assuming equilibrium and false otherwise.
    """
    inEquilibrium               ::  Bool

    """
    A boolean, which is true in case of assuming recombination processes and false otherwise.
    """
    recombinationOn             ::  Bool

    """
    A boolean, which is true in case of additionally defining conditions at the inner interfaces.
    """
    innerInterfaces             ::  Bool


    # number of boundary regions

    """
    An array for the given applied voltages at the contacts.
    """
    contactVoltage              ::  Array{Float64,1}

    """
    An array for the given Fermi level at the contacts.
    """
    bFermiLevel                 ::  Array{Float64,1}


    # number of carriers

    """
    An array with the corresponding charge numbers ``z_\\alpha`` for all carriers ``\\alpha``.
    """
    chargeNumbers               ::  Array{Float64,1}

    """
    An array with the corresponding distribution function ``\\mathcal{F}_\\alpha`` for all carriers ``\\alpha``.
    """
    F                           ::  Array{Function,1}


    # number of boundary regions x number of carriers

    """
    A 2D array with the corresponding boundary band-edge energy values ``E_\\alpha`` for each carrier ``\\alpha``.
    """
    bBandEdgeEnergy             ::  Array{Float64,2}

    """
    A 2D array with the corresponding boundary effective density of states values ``N_\\alpha`` for each carrier ``\\alpha``.
    """
    bDensityOfStates            ::  Array{Float64,2}

    """
    A 2D array with the corresponding boundary doping values for each carrier ``\\alpha``.
    """
    bDoping                     ::  Array{Float64,2}

    """
    A 2D array with the corresponding boundary velocity values for each carrier ``\\alpha``,
    when assuming Schottky contacts.
    """
    bVelocity                   ::  Array{Float64,2}


    # number of regions x number of carriers
    
    """
    A 2D array with the corresponding doping values for each carrier ``\\alpha`` on each region.
    """
    doping                      ::  Array{Float64,2}

    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha`` for each carrier ``\\alpha`` on each region.
    """
    densityOfStates             ::  Array{Float64,2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha`` for each carrier ``\\alpha`` on each region.
    """
    bandEdgeEnergy              ::  Array{Float64,2}

    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha`` for each carrier ``\\alpha`` on each region.
    """
    mobility                    ::  Array{Float64,2}


    # number of regions x 2 (for electrons and holes only!)

    """
    A 2D array with the corresponding SRH lifetimes ``\\tau_n, \\tau_p`` for electrons and holes.
    """
    recombinationSRHLifetime    ::  Array{Float64,2}

    """
    A 2D array with the corresponding SRH trap densities ``n_{\\tau}, p_{\\tau}`` for electrons and holes.
    """
    recombinationSRHTrapDensity ::  Array{Float64,2}

    """
    A 2D array with the corresponding Auger coefficients for electrons and holes.
    """
    recombinationAuger          ::  Array{Float64,2}


    # number of regions

    """
    A region dependent dielectric constant.
    """
    dielectricConstant          ::  Array{Float64,1}

    """
    A region dependent array for the emitted light in the generation process.
    """
    generationEmittedLight      ::  Array{Float64,1}

    """
    A region dependent array for the prefactor in the generation process.
    """
    generationPrefactor         ::  Array{Float64,1}

    """
    A region dependent array for the absorption coefficient in the generation process.
    """
    generationAbsorption        ::  Array{Float64,1}

    """
    A region dependent array for the radiative recombination rate.
    """
    recombinationRadiative      ::  Array{Float64,1}
    

    # number of nodes
    """
    A node dependent dielectric constant.
    """
    dielectricConstantNode      ::  Array{Float64,1}


    # number of nodes x number of carriers

    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    mobilityNode                ::  Array{Float64,2}

    """
    A 2D array with the corresponding doping values for each carrier ``\\alpha`` on each node.
    """
    dopingNode                  ::  Array{Float64,2}

    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    densityOfStatesNode         ::  Array{Float64,2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    bandEdgeEnergyNode          ::  Array{Float64,2}

    # standard constructor
    # ChargeTransportData(... all args ...) = new(... all args ...)

end



function emptyFunction()
end


"""
$(TYPEDSIGNATURES)

Simplified constructors for ChargeTransportData which only take the
number of regions, number of boundary regions and the number
of charge carriers as input.
Note that with these constructors it is also possible to simulate graded interfaces.
Almost all physical parameters can be chosen region or node dependent.
"""
function ChargeTransportData(numberOfNodes::Int64, numberOfRegions=3::Int64, numberOfBoundaryRegions=2::Int64, ;numberOfSpecies=3 ::Int64)
    # DA: Currently, we just initialize node and region dependent data with spzero and the user is setting one of them in the main file.
    # Could this be implemented better or is this idea enough?
    ChargeTransportData(

    # integer numbers
    numberOfNodes,
    numberOfRegions,
    numberOfBoundaryRegions,
    numberOfSpecies - 1,                                                      # number of carriers
    0,                                                                        # numberOfInterfaceSpecies

    # real numbers
    300 * K,                                                                  # temperature
    (kB * 300 * K ) / q,                                                      # thermal voltage
    0.0,                                                                      # reference energy
    0.27,                                                                     # parameter for Blakemore statistics
    1.0,                                                                      # λ1: embedding parameter for NLP
    0.0,                                                                      # λ2: embedding parameter for G
    1.0,                                                                      # λ3: embedding parameter for electro chemical reaction
    0.0,                                                                      # r0 prefactor from electro-chemical reaction

    # booleans
    true,                                                                     # inEquilibrium
    true,                                                                     # recombinationOn
    false,                                                                    # innerInterfaces

    # number of boundary regions
    Array{Float64,1}(undef, numberOfBoundaryRegions),                         # contactVoltage
    Array{Float64,1}(undef, numberOfBoundaryRegions),                         # Fermi level at boundary

    # number of charge carriers = number of species - 1
    Array{Float64,1}(undef, numberOfSpecies-1),                               # chargeNumbers
    fill!(similar(Array{Function,1}(undef, numberOfSpecies-1),Function),exp), # F (Boltzmann)

    # number of carriers x number of boundary regions
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bBandEdgeEnergy
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bDensityOfStates
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bDoping
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # velocity at boundary; Schottky contacts

    # number of charge carriers x number of regions
    spzeros(Float64,  numberOfSpecies-1,numberOfRegions),                     # doping  
    spzeros(Float64, numberOfSpecies-1,numberOfRegions),                      # densityOfStates
    spzeros(Float64, numberOfSpecies-1,numberOfRegions),                      # bandEdgeEnergy
    spzeros(Float64, numberOfSpecies-1,numberOfRegions),                      # mobility
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationSRHLifetime
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationSRHTrapDensity
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationAuger

    # number of regions
    spzeros(Float64, numberOfRegions),                                        # dielectricConstant
    spzeros(Float64, numberOfRegions),                                        # generationEmittedLight
    spzeros(Float64, numberOfRegions),                                        # generationPrefactor
    spzeros(Float64, numberOfRegions),                                        # generationAbsorption
    spzeros(Float64, numberOfRegions),                                        # recombinationRadiative

    # number of nodes
    spzeros(Float64, numberOfNodes),                                          # dielectricConstantNode  

    # number of carriers x number of nodes
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # mobilityNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # dopingNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # densityOfStatesNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes)                         # bandEdgeEnergyNode
    )    
    
end


function Base.show(io::IO, this::ChargeTransportData)
    for name in fieldnames(typeof(this))[1:end-5] # exclude the nodal dependent values
        @printf("%30s = ",name)
        println(io,getfield(this,name))
    end
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for interior nodes.
"""
function etaFunction(u, node::VoronoiFVM.Node, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64)
    E  = data.bandEdgeEnergy[icc, node.region] + data.bandEdgeEnergyNode[icc, node.index]
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for boundary nodes.
"""
function etaFunction(u, bnode::VoronoiFVM.BNode, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64)
    # bnode.index refers to index in overall mesh
    E  = data.bBandEdgeEnergy[icc, bnode.region] + data.bandEdgeEnergyNode[icc, bnode.index]
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for edges.
"""

function etaFunction(u, edge::VoronoiFVM.Edge, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64, nodeEdge)
    E  = data.bandEdgeEnergy[icc, edge.region] + data.bandEdgeEnergyNode[icc, nodeEdge] # if I do not put +1, I run into bounds error. It seems that VoronoiFVM allows edge.index = 0?
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
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
function breactionOhmic!(f, u, bnode, data)

    # parameters
    α          = 1.0 / VoronoiFVM.Dirichlet       # tiny penalty value
    α          = 1.0e-10                          # tiny penalty value
    ipsi       = data.numberOfCarriers + 1        # final index for electrostatic potential

    # NICHT SCHÖN: Problem interior and boundary nodes sind beide bnodes...
    if bnode.region == 1 || bnode.region == 2 

        for icc = 1:data.numberOfCarriers

            eta     = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)

            f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * (data.bDoping[icc, bnode.region] + data.dopingNode[icc, bnode.index])                             # subtract doping
            f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * (data.bDensityOfStates[icc, bnode.region] + data.densityOfStatesNode[icc, bnode.index]) * data.F[icc](eta)  # add charge carrier

            # boundary conditions for charge carriers are set in main program
            f[icc]  = 0.0

        end
    if data.innerInterfaces == true
        elseif bnode.region == 3 || bnode.region == 4 
            # NICHT SCHÖN: Problem interior and boundary nodes sind beide bnodes...  
            iphia             = 3
            iphiaj1, iphiaj2  = 5:6

            E1                = data.bBandEdgeEnergy[iphia, 3] 
            E2                = data.bBandEdgeEnergy[iphia, 4] 
            DOS1              = data.bDensityOfStates[iphia, 3] 
            DOS2              = data.bDensityOfStates[iphia, 4] 
            C01               = data.bDoping[iphia, 3]
            C02               = data.bDoping[iphia, 4]

            β                 = 0.5 # can be between 0 and 1 
            κ                 = 1 # either 0 or 1
            r0                = data.r0

            if bnode.region == 3
                etaInterfaceAnion = data.chargeNumbers[iphia] / data.UT * ( (u[iphiaj1] - u[ipsi]) + E1 / q )
    
                f[ipsi]        =  - q  *  ( data.chargeNumbers[iphia] * DOS1^(2/3) * data.F[iphia](etaInterfaceAnion) - C01^(2/3) ) # (1.4.5) @ left inner boundary 
    
                if data.inEquilibrium == true
                    f[iphia]   = u[iphia]
                    f[iphiaj1] = u[iphiaj1]
                else
    
                f[iphia]       =   data.λ3 * q * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj1, ipsi, β, κ, DOS1, E1) ) # (1.4.8) @ left inner boundary 
                f[iphiaj1]     = - data.λ3 * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj1, ipsi, β, κ, DOS1, E1) ) # (1.4.7) @ left inner boundary (right-hand side of equation)
    
                end

            elseif bnode.region == 4
                etaInterfaceAnion = data.chargeNumbers[iphia] / data.UT * ( (u[iphiaj2] - u[ipsi]) + E2 / q )
    
                f[ipsi]        =  -  q *  ( data.chargeNumbers[iphia] * DOS2^(2/3) * data.F[iphia](etaInterfaceAnion) - C02^(2/3) ) # (1.4.5) @ rigth inner boundary 
    
    
                if data.inEquilibrium == true
                    f[iphia]   = u[iphia]
                    f[iphiaj2] = u[iphiaj2]
                else
    
                f[iphia]       = - data.λ3 *  q * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj2, ipsi, β, κ, DOS2, E2) ) # (1.4.8) @ right inner boundary 
                f[iphiaj2]     = - data.λ3 * ( r0 * electrochemicalReaction(data, u, iphia, ipsi, iphiaj2, ipsi, β, κ, DOS2, E2) ) # (1.4.7) @ right inner boundary (right-hand side of equation)
                end
            end

        end # if clause for inner interface       

    end

    f[ipsi] = -1 / α *  q * data.λ1 * f[ipsi]

end

function electrochemicalReaction(data, u, iphia, ipsi, iphiaJunction, ipsiJunction, β, κ, DOS, E) # (1.4.9)
 

    etaExp             = data.chargeNumbers[iphia] / data.UT * ( (u[iphia] - u[iphiaJunction]) + E / q ) 
    expTerm            =  exp( β * etaExp ) - exp( (β - 1) * etaExp)

    etaInterfaceAnion  = data.chargeNumbers[iphia] / data.UT * ( (u[iphiaJunction] - u[ipsiJunction]) + E / q )
    etaAnion           = data.chargeNumbers[iphia] / data.UT * ( (u[iphia] - u[ipsi]) + E / q )

    densFactor         = ( (DOS^(2/3) * data.F[iphia](etaInterfaceAnion) )^(1/2) * (DOS * data.F[iphia](etaAnion) )^(- 1/2) )^κ

    return densFactor * expTerm

end

function bstorage!(f, u, bnode, data)

    iphia             = 3
    ipsi              = data.numberOfCarriers + 1        # final index for electrostatic potential
    iphiaj1, iphiaj2  = 5:6

    E1                = data.bBandEdgeEnergy[iphia, 3] 
    E2                = data.bBandEdgeEnergy[iphia, 4] 
    DOS1              = data.bDensityOfStates[iphia, 3] 
    DOS2              = data.bDensityOfStates[iphia, 4] 

    if bnode.region == 3

        # (1.4.7) @ left inner boundary (left-hand side of equation)
        etaInterfaceAnion = data.chargeNumbers[iphia] / data.UT * ( (u[iphiaj1] - u[ipsi]) + E1 / q ) 
        f[iphiaj1]        = data.chargeNumbers[iphia] * DOS1^(2/3) * data.F[iphia](etaInterfaceAnion) 

    elseif bnode.region == 4

        # (1.4.7) @ right inner boundary (left-hand side of equation)
        etaInterfaceAnion = data.chargeNumbers[iphia] / data.UT * ( (u[iphiaj2] - u[ipsi]) + E2 / q )
        f[iphiaj2]        = data.chargeNumbers[iphia] *  DOS2^(2/3) * data.F[iphia](etaInterfaceAnion)

    end
end

"""
$(TYPEDSIGNATURES)

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

[Note that this way of implementation is not well tested yet.]
"""
function breactionSchottky!(f, u, bnode, data)
    # DA: Still need to be well tested!!

    ipsi          = data.numberOfCarriers + 1        # final index for electrostatic potential

    # bnode.coord
    # if bnode.region == 1
    #     f[ipsi] = -u[ipsi] 
    # elseif bnode.region == 2
    #     f[ipsi] =  + u[ipsi] - data.λ1 *((data.bFermiLevel[2] - data.bFermiLevel[1])/q ) - data.contactVoltage[bnode.region]  
    # end

    for icc = 1:data.numberOfCarriers-1
       
        if bnode.region == 1
            phi = data.bFermiLevel[1]
        elseif bnode.region == 2
            phi = data.bFermiLevel[2]
        end
        if bnode.region == 1 || bnode.region == 2
            E          = data.bBandEdgeEnergy[icc, bnode.region] + data.bandEdgeEnergyNode[icc, bnode.index]
        etaFix = data.chargeNumbers[icc] / data.UT * (  (-phi + E ) / q  )
        eta    = data.chargeNumbers[icc] / data.UT * (  (u[icc]  - u[ipsi]) + E / q )

        f[icc] =  data.chargeNumbers[icc] * q * data.λ1* data.bVelocity[icc, bnode.region] * (  (data.bDensityOfStates[icc, bnode.region] + data.densityOfStatesNode[icc, bnode.index])  * (data.F[icc](eta) - data.F[icc](etaFix)  ))
        end

        # if icc == 3
        #     if bnode.region == 4
        #         E          = data.bBandEdgeEnergy[icc, bnode.region] + data.bandEdgeEnergyNode[icc, bnode.index]
        #         phi = phi_right
        #     etaFix = data.chargeNumbers[icc] / data.UT * (  (-phi + E ) / q  )
        #     eta    = data.chargeNumbers[icc] / data.UT * (  (u[icc]  - u[ipsi]) + E / q )
    
        #     f[icc] =  data.chargeNumbers[icc] * q * data.λ1* 0.0 * (  (data.bDensityOfStates[icc, bnode.region] + data.densityOfStatesNode[icc, bnode.index])  * (data.F[icc](eta) - data.F[icc](etaFix)  ))
        #     end
    
        # end

    end

end

"""
$(TYPEDSIGNATURES)

The generation rate ``G``, which occurs in the right-hand side of the
continuity equations. Currently, we assume a region dependent constant value.
"""
function generation(data, ireg)

    return data.λ2 * data.generationEmittedLight[ireg]    # Phil considers a uniform generation rate (but only in the intrinsic layer)

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
function reaction!(f, u, node, data)

    # indices
    iphin                 = 1
    iphip                 = 2
    ipsi                  = data.numberOfCarriers + 1            # final index for electrostatic potential
    ireg                  = node.region
    inode                 = node.index

    # rhs of NLP (charge density)
    for icc = 1:data.numberOfCarriers
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * (data.doping[icc, node.region] + data.dopingNode[icc, node.index])  # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * (data.densityOfStates[icc, node.region] + data.densityOfStatesNode[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    # rhs of continuity equations for electron and holes (bipolar reaction)
    
    if data.inEquilibrium == true 

        for icc = 1:data.numberOfCarriers
         f[icc] = u[icc] - 0.0
        end

    else
        if data.recombinationOn == true  
            n                     = computeDensities(u, data, inode, ireg, iphin, ipsi, true)  # true for interior region
            p                     = computeDensities(u, data, inode, ireg, iphip, ipsi, true) 
            exponentialTerm       = exp((q *u[iphin] - q  * u[iphip]) / (kB*data.temperature))
            excessCarrierDensTerm = n*p * (1.0 - exponentialTerm)

            for icc in [iphin, iphip] 

             
                # radiative recombination
                kernelRadiative = data.recombinationRadiative[ireg]
                
                # Auger recombination
                kernelAuger     = (data.recombinationAuger[iphin, ireg] * n + data.recombinationAuger[iphip, ireg] * p)
                
                # SRH recombination
                kernelSRH       = 1.0 / (  data.recombinationSRHLifetime[iphip, ireg] * (n + data.recombinationSRHTrapDensity[iphin, ireg]) + data.recombinationSRHLifetime[iphin, ireg] * (p + data.recombinationSRHTrapDensity[iphip, ireg]) )
                
                # full recombination
                f[icc]          = q* data.chargeNumbers[icc]* (kernelRadiative + kernelAuger + kernelSRH)*  excessCarrierDensTerm  - q * data.chargeNumbers[icc] * generation(data, ireg)
            end

        else
            for icc = 1: data.numberOfCarriers
                f[icc]          = 0.0 #- q * data.chargeNumbers[icc] * generation(data, ireg) 
            end
        end

        for icc in iphip+1:data.numberOfCarriers
            f[icc]              = u[icc] - 0.0
        end
    
    end
    
    f[ipsi]                     = - q * data.λ1 * f[ipsi]
end


"""
$(TYPEDSIGNATURES)

The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have 

``f[n_\\alpha] =  z_\\alpha  q ∂_t n_\\alpha`` 

and for the electrostatic potential
``f[ψ] = 0``.
[Note that this method is not tested yet!]
"""
function storage!(f, u, node, data)
    ipsi       = data.numberOfCarriers + 1
    
    for icc = 1:data.numberOfCarriers

        eta    = etaFunction(u, node, data, icc, ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = q * data.chargeNumbers[icc] * (data.densityOfStates[icc, node.region] + data.densityOfStatesNode[icc, node.index]) * data.F[icc](eta)
    end

    f[ipsi] = 0.0
end

"""
$(SIGNATURES)

Compute trap densities for a given trap energy.
[Currently, only done for the Boltzmann statistics.]

"""
function trapDensity(icc, ireg, data, Et) 
    data.densityOfStates[icc, ireg] * exp( data.chargeNumbers[icc] * (data.bandEdgeEnergy[icc, ireg] - Et) / (kB * data.temperature)) # need to subtract Eref
end


"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme.

"""
function ScharfetterGummel!(f, u, edge, data)
    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    nodel    = edge.node[2]
    nodek    = edge.node[1]
    
    dpsi     = ul[ipsi] - uk[ipsi]
    f[ipsi]  = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge)
        etal               = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi, nodeEdge)

        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        bp, bm             = fbernoulli_pm(data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT)
        f[icc]             = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
    
    end

end

"""
$(TYPEDSIGNATURES)

The classical Scharfetter-Gummel flux scheme for 
possible space-dependent DOS and band-edge energies. For these parameters
the discretization scheme is modified. [insert continuous flux etc ...]

"""
function ScharfetterGummelGraded!(f, u, edge, data)
    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    nodel    = edge.node[2]
    nodek    = edge.node[1]
    
    dpsi     = ul[ipsi] - uk[ipsi]
    dpsiEps  = (data.dielectricConstant[ireg]  + (data.dielectricConstantNode[nodel] + data.dielectricConstantNode[nodek])/2) * dpsi
    f[ipsi]  = - ε0 * dpsiEps
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.UT

        etak               = data.chargeNumbers[icc] / data.UT * ( (uk[icc] - uk[ipsi]) + (data.bandEdgeEnergyNode[icc, nodek] + data.bandEdgeEnergy[icc, ireg]) / q )
        etal               = data.chargeNumbers[icc] / data.UT * ( (ul[icc] - ul[ipsi]) + (data.bandEdgeEnergyNode[icc, nodel] + data.bandEdgeEnergy[icc, ireg]) / q )


        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]
        mobility           = data.mobility[icc, ireg] + (data.mobilityNode[icc, nodel] + data.mobilityNode[icc, nodek])/2
        densityOfStatesl   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodel])
        densityOfStatesk   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodek])
        
        if data.densityOfStatesNode[icc, nodel] ≈ 0.0 || data.densityOfStatesNode[icc, nodek] ≈ 0.0
            bp, bm         = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT ) 
        else
            bp, bm         = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT - (log(data.densityOfStatesNode[icc, nodel]) -log(data.densityOfStatesNode[icc,nodek])) ) 
        end

        f[icc]             =  -j0  * mobility * ( bm  * densityOfStatesl * data.F[icc](etal) - bp *  densityOfStatesk * data.F[icc](etak) )

    end

end


"""
$(TYPEDSIGNATURES)

The Sedan flux scheme.

"""
function Sedan!(f, u, edge, data)
    
    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    nodel    = edge.node[2]
    nodek    = edge.node[1]
    
    dpsi     = ul[ipsi] - uk[ipsi]
    f[ipsi]  = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u, edge, data, icc, ipsi, nodeEdge)
        etal               = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u, edge, data, icc, ipsi, nodeEdge)

        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        Q                  = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
end

end


"""
$(TYPEDSIGNATURES)

The Sedan flux scheme for 
possible space-dependent DOS and band-edge energies. For these parameters
the discretization scheme is modified. [insert continuous flux etc ...]

"""
function SedanGraded!(f, u, edge, data)

    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    nodel    = edge.node[2]
    nodek    = edge.node[1]
    
    dpsi     = ul[ipsi] - uk[ipsi]
    dpsiEps  = (data.dielectricConstant[ireg]  + data.dielectricConstantNode[nodel]) * ul[ipsi] - (data.dielectricConstant[ireg] + data.dielectricConstantNode[nodek]) * uk[ipsi]
    f[ipsi]  = - ε0 * dpsiEps
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 = data.chargeNumbers[icc] * q * data.UT
        etak               = data.chargeNumbers[icc] / data.UT * ( (uk[icc] - uk[ipsi]) + (data.bandEdgeEnergyNode[icc, nodek] + data.bandEdgeEnergy[icc, ireg]) / q )
        etal               = data.chargeNumbers[icc] / data.UT * ( (ul[icc] - ul[ipsi]) + (data.bandEdgeEnergyNode[icc, nodel] + data.bandEdgeEnergy[icc, ireg]) / q )

        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]
        mobilityl          = (data.mobility[icc, ireg] + data.mobilityNode[icc, nodel])
        mobilityk          = (data.mobility[icc, ireg] + data.mobilityNode[icc, nodek])
        densityOfStatesl   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodel])
        densityOfStatesk   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodek])

        if data.densityOfStatesNode[icc, nodel] ≈ 0.0 || data.densityOfStatesNode[icc, nodek] ≈ 0.0
            Q              = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )

        else
            Q              = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) - (log(data.densityOfStatesNode[icc, nodel]) -log(data.densityOfStatesNode[icc,nodek])) )

        end

        bp, bm         = fbernoulli_pm(Q)
        f[icc]         = - j0  * ( bm  * mobilityl * densityOfStatesl * data.F[icc](etal) - bp * mobilityk * densityOfStatesk * data.F[icc](etak) )
end

end


"""
$(TYPEDSIGNATURES)

The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is 
used for the regularization of the removable singularity.

"""
function diffusionEnhanced!(f, u, edge, data)
    tolReg  = 1.0e-13;
    
    uk      = viewK(edge, u)
    ul      = viewL(edge, u)
    
    ipsi    = data.numberOfCarriers + 1
    ireg    = edge.region
    nodel   = edge.node[2]
    nodek   = edge.node[1]
    
    dpsi    = ul[ipsi] - uk[ipsi]
    f[ipsi] = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

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

        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        bp, bm             = fbernoulli_pm(data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / (data.UT * g))
        f[icc]             = - j0 * g * (  bm * data.F[icc](etal) - bp * data.F[icc](etak))
    end

end

"""
$(TYPEDSIGNATURES)

The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation which arise
when considering the generalized Scharfetter-Gummel scheme in case of Blakemore statistics.
Hence, it should be exclusively worked with, when considering the Blakemore distribution.

"""
function KopruckiGaertner!(f, u, edge, data)
    max_iteration = 200          # for Newton solver
    it            = 0            # number of iterations (newton)
    damp          = 0.1          # damping factor
    
    uk            = viewK(edge, u)
    ul            = viewL(edge, u)
    nodel         = edge.node[2]
    nodek         = edge.node[1]
    
    ipsi          = data.numberOfCarriers + 1
    ireg          = edge.region
    
    dpsi          = ul[ipsi] - uk[ipsi]
    f[ipsi]        =  - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers 

        j0                  = data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak                = etaFunction(uk, edge, data, icc, ipsi, nodek) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)
        etal                = etaFunction(ul, edge, data, icc, ipsi, nodel) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi, nodeEdge)

        bandEdgeDifference  = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        # use Sedan flux as starting guess
        Q                   = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm              = fbernoulli_pm(Q)
        jInitial            = ( bm * data.F[icc](etal)  - bp * data.F[icc](etak))

        implicitEq(j::Real) = (fbernoulli_pm(data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q)) /data.UT + data.γ * j )[2] * exp(etal) - fbernoulli_pm(data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q) /data.UT) - data.γ * j)[1] * exp(etak)) - j

        delta               = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
        oldup = 1.0
        while (it < max_iteration)
            Fval     = implicitEq(jInitial)
            dFval    = ForwardDiff.derivative(implicitEq, jInitial)

            if isnan(value(dFval)) || value(abs(dFval)) < delta
                @show value(jInitial), value(Fval), value(dFval)
                error("singular derivative")
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
end


"""
$(TYPEDSIGNATURES)

The argument of the distribution function

``z_\\alpha / U_T  ( (\\varphi_\\alpha - \\psi) + E_\\alpha / q )``

for floats.
"""
function etaFunction(u, data::VoronoiFVM.AbstractData, node, region, icc::Int64, ipsi::Int64, in_region::Bool)

    if in_region == true
        E  = data.bandEdgeEnergy[icc, region] + data.bandEdgeEnergyNode[icc, node]
    elseif in_region == false
        E  = data.bBandEdgeEnergy[icc, region] + data.bandEdgeEnergyNode[icc, node]
    end

    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities.

"""
function computeDensities(u, data::VoronoiFVM.AbstractData, node, region, icc::Int, ipsi::Int, in_region::Bool)

    if in_region == false
        (data.bDensityOfStates[icc, region] + data.densityOfStatesNode[icc, node] ) * data.F[icc](etaFunction(u, data, node, region, icc, ipsi, in_region::Bool))
    elseif in_region == true
        (data.densityOfStates[icc, region] + data.densityOfStatesNode[icc, node])* data.F[icc](etaFunction(u, data, node, region, icc, ipsi, in_region::Bool)) 
    end
        
end


"""

$(TYPEDSIGNATURES)

For given potentials in vector form, compute corresponding vectorized densities.
[Caution: this was not tested for multidimensions.]
"""
function computeDensities(grid, data, sol)
    ipsi         = data.numberOfCarriers + 1 
    densities    = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
        
    bfacenodes   = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    cellRegions  = copy(grid[CellRegions])
    cellRegions  = push!(cellRegions, grid[CellRegions][end]) #  enlarge region by final cell
    
    if dim_space(grid) > 1
        println("ComputeDensities is so far only tested in 1D")
    end
    
    for icc in 1:data.numberOfCarriers

        for node in 1:data.numberOfNodes
            in_region = true
            u         = sol[:, node]
            region    = cellRegions[node]

            if node in bfacenodes
                in_region = false
                indexNode = findall(x -> x == node, vec(bfacenodes))[1]  # we need to know which index the node has in bfacenodes
                region    = bfaceregions[indexNode]                      # since the corresponding region number is at the same index
            end

            densities[icc, node] = computeDensities(u, data, node, region, icc, ipsi, in_region)
        end
    
    end
    
    return densities

end


"""

$(SIGNATURES)

For given solution in vector form, compute corresponding vectorized band-edge energies and Fermi level.
[Caution: this was not tested for multidimensions.]
"""
function computeEnergies(grid, data, sol)
    
    ipsi         = data.numberOfCarriers + 1
    energies     = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
    fermiLevel   = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
    
    cellregions  = grid[CellRegions]
    cellregions  = push!(cellregions, cellregions[end])
    
    for icc in 1:data.numberOfCarriers

        for inode in 1:data.numberOfNodes
             E                      = data.bandEdgeEnergy[icc, cellregions[inode]] + data.bandEdgeEnergyNode[icc, inode]
             energies[icc, inode]   = E - q * sol[ipsi, inode]
             fermiLevel[icc, inode] = -q * sol[icc, inode]
        end
    
    end
    
    return energies, fermiLevel

end


"""

$(TYPEDSIGNATURES)

Compute the electro-neutral solution for the Boltzmann approximation. 
It is obtained by setting the left-hand side in
the Poisson equation equal to zero and solving for ``\\psi``.
The charge carriers may obey different statitics functions.
Currently, this one is not well tested for the case of charge carriers beyond electrons and holes.
"""
function electroNeutralSolution!(data, grid; Newton=false)

    if data.numberOfCarriers > 2
        error("this method is currently only working for electrons and holes")
    end

    solution        = zeros(length(grid[Coordinates]))
    iccVector       = collect(1:data.numberOfCarriers)
    zVector         = data.chargeNumbers[iccVector]
    FVector         = data.F[iccVector]
    regionsAllCells = copy(grid[CellRegions])
    regionsAllCells = push!(regionsAllCells, grid[CellRegions][end]) #  enlarge region by final cell
    phi             = 0.0                                            # in equilibrium set to 0
    psi0_initial    = 0.5

    for index = 1:length(regionsAllCells) - 1
        
        ireg          = regionsAllCells[index]
        zVector       = data.chargeNumbers[iccVector]
        FVector       = data.F[iccVector]
        regionsOfCell = regionsAllCells[grid[CellNodes][:,index]]   # all regions of nodes belonging to cell for given index

        # average following quantities if needed among all regions
        EVector = Float64[]; CVector = Float64[]; NVector = Float64[]

        for icc = 1:data.numberOfCarriers
            push!(EVector, sum(data.bandEdgeEnergy[icc, regionsOfCell])  / length(regionsOfCell) + data.bandEdgeEnergyNode[icc,index])
            push!(CVector, sum(data.doping[icc, regionsOfCell])          / length(regionsOfCell) + data.dopingNode[icc, index])
            push!(NVector, sum(data.densityOfStates[icc, regionsOfCell]) / length(regionsOfCell) + data.densityOfStatesNode[icc, index])
        end
        # rhs of Poisson's equation as anonymous function depending on psi0
        f = psi0 -> chargeDensity(psi0, phi, data.UT, EVector, zVector, CVector, NVector, FVector)

        if !Newton
            try
                solution[index + 1] = fzero(f, psi0_initial)
            catch 
                psi0_initial        = 2.0
                solution[index + 1] = fzero(f, psi0_initial)
                psi0_initial        = 0.25
            end
        else 
            D(f)                    = psi0 -> ForwardDiff.derivative(f, float(psi0))
            solution[index + 1]     = find_zero((f, D(f)), psi0_initial)
        end
    end

    # fill in last values, same as second to last
    solution[1] = solution[2]

    return solution

end


"""

$(TYPEDSIGNATURES)

Compute the charge density, i.e. the right-hand side of Poisson's equation.

"""
function chargeDensity(psi0, phi, UT, EVector, chargeNumbers, dopingVector, dosVector, FVector)
    # https://stackoverflow.com/questions/45667291/how-to-apply-one-argument-to-arrayfunction-1-element-wise-smartly-in-julia
    sum(-chargeNumbers .* dopingVector) + sum(chargeNumbers .* dosVector .* (etaFunction(psi0, phi, UT, EVector, chargeNumbers) .|> FVector))
end

"""

$(TYPEDSIGNATURES)

First try of debugger. Print the Jacobi matrix for a given node, i.e.
the number of the node in the grid and not the excact coordinate.
This is only done for the one dimensional case so far.
[insert some information about VoronoiFVM]
"""
function printJacobi(node, sys)
    ct_data = data(sys)
    numberOfNodes = ct_data.numberOfNodes
    if node == 1
        println(sys.matrix[1:3, 1:9])
    elseif node == numberOfNodes
        println(sys.matrix[3*numberOfNodes-2:3*numberOfNodes, 3*numberOfNodes-8:3*numberOfNodes])
    else
        println(sys.matrix[3*node-2:3*node, 3*node-5:3*node+3])
    end
end
