"""
$(TYPEDEF)

Struct holding physical data for drift-diffusion simulation of semiconductor device.
If there are ``N`` number of species, it is assumed that the first ``N-1``ones
correspond to the charge carriers and the final one to the electrostatic potential.

$(TYPEDFIELDS)
"""

mutable struct ChargeTransportData <: VoronoiFVM.AbstractData

    # integer numbers
    numberOfNodes               ::  Int64
    numberOfRegions             ::  Int64
    numberOfBoundaryRegions     ::  Int64
    numberOfCarriers            ::  Int64

    # real numbers
    temperature                 ::  Float64
    UT                          ::  Float64
    Eref                        ::  Float64
    γ                           ::  Float64
    λ1                          ::  Float64
    λ2                          ::  Float64

    # booleans
    inEquilibrium               ::  Bool
    recombinationOn             ::  Bool

    # number of boundary regions
    contactVoltage              ::  Array{Float64,1}

    # number of carriers
    chargeNumbers               ::  Array{Float64,1}
    F                           ::  Array{Function,1}

    # number of boundary regions x number of carriers
    bBandEdgeEnergy             ::  Array{Float64,2}
    bDensityOfStates            ::  Array{Float64,2}
    bDoping                     ::  Array{Float64,2}

    # number of regions x number of carriers
    doping                      ::  Array{Float64,2}
    densityOfStates             ::  Array{Float64,2}
    bandEdgeEnergy              ::  Array{Float64,2}
    mobility                    ::  Array{Float64,2}
    recombinationSRHLifetime    ::  Array{Float64,2}
    recombinationSRHTrapDensity ::  Array{Float64,2}
    recombinationAuger          ::  Array{Float64,2}

    # number of regions
    dielectricConstant          ::  Array{Float64,1}
    generationEmittedLight      ::  Array{Float64,1}
    generationPrefactor         ::  Array{Float64,1}
    generationAbsorption        ::  Array{Float64,1}
    recombinationRadiative      ::  Array{Float64,1}
    
    # number of nodes
    dielectricConstantNode      ::  Array{Float64,1}

    # number of nodes x number of carriers
    mobilityNode                ::  Array{Float64,2}
    dopingNode                  ::  Array{Float64,2}
    densityOfStatesNode         ::  Array{Float64,2}   # still needs to be implemented
    bandEdgeEnergyNode          ::  Array{Float64,2}   # still needs to be implemented

    # standard constructor
    # ChargeTransportData(... all args ...) = new(... all args ...)

end

function emptyFunction()
end

"""

$(SIGNATURES)

Simplified constructors for ChargeTransportData which takes only the
number of regions, number of boundary regions and the number
of charge carriers as input.

"""

function ChargeTransportData(numberOfNodes::Int64, numberOfRegions=3::Int64, numberOfBoundaryRegions=2::Int64, ;numberOfSpecies=3 ::Int64)
    ChargeTransportData(

    # integer numbers
    numberOfNodes,
    numberOfRegions,
    numberOfBoundaryRegions,
    numberOfSpecies - 1,                                                     # number of carriers

    # real numbers
    300 * K,                                                                 # temperature
    (kB * 300 * K ) / q,                                                     # thermal voltage
    0.0,                                                                     # reference energy
    0.27,                                                                    # parameter for Blakemore statistics
    1.0,                                                                     # λ1: embedding parameter for NLP
    0.0,                                                                     # λ2: embedding parameter for G

    # booleans
    true,                                                                    # inEquilibrium
    true,                                                                    # recombinationOn

    # number of boundary regions
    Array{Float64,1}(undef, numberOfBoundaryRegions),                        # contactVoltage

    # number of charge carriers = number of species - 1
    Array{Float64,1}(undef, numberOfSpecies-1),                               # chargeNumbers
    fill!(similar(Array{Function,1}(undef, numberOfSpecies-1),Function),exp), # F (Boltzmann)

    # number of carriers x number of boundary regions
    Array{Float64,2}(undef, numberOfSpecies-1, numberOfBoundaryRegions),      # bBandEdgeEnergy
    Array{Float64,2}(undef, numberOfSpecies-1, numberOfBoundaryRegions),      # bDensityOfStates
    zeros(Float64,          numberOfSpecies-1, numberOfBoundaryRegions),      # bDoping

    # number of charge carriers x number of regions
    zeros(Float64,      numberOfSpecies-1,numberOfRegions),                   # doping

    Array{Float64,2}(undef, numberOfSpecies-1,numberOfRegions),               # densityOfStates
    Array{Float64,2}(undef, numberOfSpecies-1,numberOfRegions),               # bandEdgeEnergy
    Array{Float64,2}(undef, numberOfSpecies-1,numberOfRegions),               # mobility
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationSRHLifetime
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationSRHTrapDensity
    Array{Float64,2}(undef,2, numberOfRegions),                               # recombinationAuger

    # number of regions
    Array{Float64,1}(undef ,numberOfRegions),                                 # dielectricConstant
    Array{Float64,1}(undef, numberOfRegions),                                 # generationEmittedLight
    Array{Float64,1}(undef, numberOfRegions),                                 # generationPrefactor
    Array{Float64,1}(undef, numberOfRegions),                                 # generationAbsorption
    Array{Float64,1}(undef, numberOfRegions),                                 # recombinationRadiative

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

$(SIGNATURES)

Simplified constructors for ChargeTransportDataGraded which takes only the
number of nodes, number of boundary regions and the number
of charge carriers as input. This is for graded structs, i.e. for devices.
where graded junctions are taken into account.

"""

mutable struct ChargeTransportDataGraded <: VoronoiFVM.AbstractData

    # integer numbers
    numberOfNodes               ::  Int64
    numberOfRegions             ::  Int64
    numberOfBoundaryRegions     ::  Int64
    numberOfCarriers            ::  Int64

    # real numbers
    temperature                 ::  Float64
    UT                          ::  Float64
    Eref                        ::  Float64
    γ                           ::  Float64
    λ1                          ::  Float64
    λ2                          ::  Float64

    # booleans
    inEquilibrium               ::  Bool
    recombinationOn             ::  Bool

    # number of boundary regions
    contactVoltage              ::  Array{Float64,1}

    # number of carriers
    chargeNumbers               ::  Array{Float64,1}
    F                           ::  Array{Function,1}

    # number of boundary regions x number of carriers
    bBandEdgeEnergy             ::  Array{Float64,2}
    bDensityOfStates            ::  Array{Float64,2}
    bDoping                     ::  Array{Float64,2}

    # number of nodes x number of carriers
    doping                      ::  Array{Float64,2}
    densityOfStates             ::  Array{Float64,2}
    bandEdgeEnergy              ::  Array{Float64,2}
    mobility                    ::  Array{Float64,2}
    recombinationSRHLifetime    ::  Array{Float64,2}
    recombinationSRHTrapDensity ::  Array{Float64,2}
    recombinationAuger          ::  Array{Float64,2}

    # number of regions
    dielectricConstant          ::  Array{Float64,1}
    generationEmittedLight      ::  Array{Float64,1}
    generationPrefactor         ::  Array{Float64,1}
    generationAbsorption        ::  Array{Float64,1}
    recombinationRadiative      ::  Array{Float64,1}
    # number of nodes
    dielectricConstantNode      ::  Array{Float64,1}
    recombinationRadiativeNode      ::  Array{Float64,1}

    # number of nodes x number of carriers
    mobilityNode                ::  Array{Float64,2}
    dopingNode                  ::  Array{Float64,2}
    densityOfStatesNode         ::  Array{Float64,2}   # still needs to be implemented
    bandEdgeEnergyNode          ::  Array{Float64,2}   # still needs to be implemented

    # standard constructor
    # ChargeTransportData(... all args ...) = new(... all args ...)

end


function ChargeTransportDataGraded(numberOfNodes::Int64, numberOfRegions=5::Int64, numberOfBoundaryRegions=2::Int64, ;numberOfSpecies=3 ::Int64)
    ChargeTransportDataGraded(
    # Note that all region dependent data is initialized with zero.

    # integer numbers
    numberOfNodes,
    numberOfRegions,
    numberOfBoundaryRegions,
    numberOfSpecies - 1,                                                      # number of carriers

    # real numbers
    300 * K,                                                                  # temperature
    (kB * 300 * K ) / q,                                                      # thermal voltage
    0.0,                                                                      # reference energy
    0.27,                                                                     # parameter for Blakemore statistics
    1.0,                                                                      # λ1: embedding parameter for NLP
    0.0,                                                                      # λ2: embedding parameter for G

    # booleans
    true,                                                                     # inEquilibrium
    true,                                                                     # recombinationOn

    # number of boundary regions
    Array{Float64,1}(undef, numberOfBoundaryRegions),                         # contactVoltage

    # number of charge carriers = number of species - 1
    Array{Float64,1}(undef, numberOfSpecies-1),                               # chargeNumbers
    fill!(similar(Array{Function,1}(undef, numberOfSpecies-1),Function),exp), # F (Boltzmann)

    # number of carriers x number of boundary regions
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bBandEdgeEnergy
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bDensityOfStates
    spzeros(Float64,  numberOfSpecies-1, numberOfBoundaryRegions),            # bDoping

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
    spzeros(Float64, numberOfNodes),                                          # recombinationRadiativeNode  

    # number of carriers x number of nodes
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # mobilityNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # dopingNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes),                        # densityOfStatesNode
    spzeros(Float64, numberOfSpecies-1,numberOfNodes)                         # bandEdgeEnergyNode
    )

end


"""

$(SIGNATURES)

The argument of the distribution function for interior nodes:

    z / UT  * ( (phi - psi) + E / q ).

"""
function etaFunction(u, node::VoronoiFVM.Node, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64)
    E  = data.bandEdgeEnergy[icc, node.region] + data.bandEdgeEnergyNode[icc, node.index]
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for boundary nodes:
    z / UT  * ( (phi_at_boundary - psi) + E / q ).
"""

function etaFunction(u, bnode::VoronoiFVM.BNode, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64)
    # bnode.index refers to index in overall mesh
    E  = data.bBandEdgeEnergy[icc, bnode.region] + data.bandEdgeEnergyNode[icc, bnode.index]
    data.chargeNumbers[icc] / data.UT * ( (data.contactVoltage[bnode.region] - u[ipsi]) + E / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for edges:

    z / UT  * ( (phi_at_edge - psi) + E / q ).

"""

function etaFunction(u, edge::VoronoiFVM.Edge, data::VoronoiFVM.AbstractData, icc::Int64, ipsi::Int64)
    E  = data.bandEdgeEnergy[icc, edge.region] + data.bandEdgeEnergyNode[icc, edge.index+1] # icell: Number of discretization cell the edge is invoked from
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for given psi and phi:

    z / UT  * ( (phi_at_edge - psi) + E / q ).

E and z are given as vectors. This function may be used to compute 
the charge density.    
"""

function etaFunction(psi, phi, UT, E::Array, z::Array)
    z ./ UT .* ( (phi - psi) .+ E / q )
end

"""
$(SIGNATURES)

Creates the boundary conditions via a penalty approach with penalty parameter 1/α.
For example, the right-hand side for the electrostatic potential is implemented as

    f[ipsi]  = -1/α *  q * ( (p - N_a) - (n - N_d) ),

assuming a bipolar semiconductor. In general, for some charge number `z_i`

    f[ipsi] =  -1/α *  q * sum_i { z_i * (c_i - N_i) }.

The boundary conditions for the charge carrier are set in the main file. Hence,

    f[icc] = 0

for all charge carriers `icc`.

"""
function breaction!(f, u, bnode, data)

    # parameters
    α          = 1.0 / VoronoiFVM.Dirichlet       # tiny penalty value
    α          = 1.0e-10                          # tiny penalty value
    ipsi       = data.numberOfCarriers + 1        # final index for electrostatic potential

    phi_right  = -4.15 * eV 
    phi_left   = -5.2 * eV 
    s_left     = [1.0e8*cm/s 1.0e8*cm/s]
    s_right    = [1.0e8*cm/s 1.0e8*cm/s]

    # NICHT SCHÖN: Problem interior and boundary nodes sind beide bnodes...
    ### TESTEN AUF ÄU?EREN RAND 
    # bnode.coord
    if bnode.region == 1 || bnode.region == 2 

        for icc = 1:data.numberOfCarriers

            eta     = etaFunction(u, bnode, data, icc, ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)

            f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * (data.bDoping[icc, bnode.region] + data.dopingNode[icc, bnode.index])                             # subtract doping
            f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * (data.bDensityOfStates[icc, bnode.region] + data.densityOfStatesNode[icc, bnode.index]) * data.F[icc](eta)  # add charge carrier

            # boundary conditions for charge carriers are set in main program
            f[icc]  = 0.0
            

        end
     
    # NICHT SCHÖN: Problem interior and boundary nodes sind beide bnodes...  
    ### TESTEN AUF INNEREN RAND  
    # elseif bnode.region == 3 || bnode.region == 4 

    #     phi_right  = -4.15 * eV 
    #     phi_left   = -5.2 * eV 
    #     s_left     = [1.0e8*cm/s 1.0e8*cm/s]
    #     s_right    = [1.0e8*cm/s 1.0e8*cm/s]

    #     for icc = 1:data.numberOfCarriers
    #         E      = data.bBandEdgeEnergy[bnode.region,icc] + data.bandEdgeEnergyNode[bnode.index,icc]
            
    #         if bnode.region == 3
    #             etaFix = data.chargeNumbers[icc] / data.UT * ((phi_left - u[ipsi]) + E / q )
    #             eta    = data.chargeNumbers[icc] / data.UT * (  (u[icc] - u[ipsi]) + E / q )
            
    #             f[icc] = -data.λ1* s_left[icc] * (  data.bDensityOfStates[bnode.region, icc] * (data.F[icc](eta) - data.F[icc](etaFix)  ))

    #         elseif bnode.region == 4
    #             etaFix = data.chargeNumbers[icc] / data.UT * ( (phi_right - u[ipsi]) + E / q )
    #             eta    = data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
    
    #             f[icc] = -data.λ1* s_right[icc] * (  data.bDensityOfStates[bnode.region, icc] * (data.F[icc](eta) - data.F[icc](etaFix)  ))
    #         end
    #     end
    end

    f[ipsi] = -1 / α *  q * data.λ1 * f[ipsi]


end


"""
(SIGNATURES)

Generation rate.
"""

function generation(data, ireg)

    return data.λ2 * data.generationEmittedLight[ireg]    # Phil considers a uniform generation rate (but only in the intrinsic layer)

end

"""
$(SIGNATURES)

Sets up the right-hand sides. Assuming a bipolar semiconductor
the right-hand side for the electrostatic potential becomes

  ``f[ψ]  = - q ((p - N_a) - (n - N_d) ) = - q  \\sum  z_i  (c_i - N_i) ``

and the right-hand sides for the charge carriers yields

``f[c_i] =  - z_i  q (G -  R) ``

for a charge number ``z_i`` and all charge carriers ``c_i``.
The recombination includes radiative, Auger and Shockley-Read-Hall
recombination.

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
    n                     = computeDensities(u, data, inode, ireg, iphin, ipsi, true)  # true for interior region
    p                     = computeDensities(u, data, inode, ireg, iphip, ipsi, true) 
    exponentialTerm       = exp((q *u[iphin] - q  * u[iphip]) / (kB*data.temperature))
    excessCarrierDensTerm = n*p * (1.0 - exponentialTerm)

    # rhs of NLP (charge density)
    for icc = 1:data.numberOfCarriers
        
        eta     = etaFunction(u, node, data, icc, ipsi) 
        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * (data.doping[icc, node.region] + data.dopingNode[icc, node.index])  # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.densityOfStates[icc, node.region] * data.F[icc](eta)   # add charge carrier

    end

    # rhs of continuity equations for electron and holes (bipolar reaction)
    
    if data.inEquilibrium == true 

        for icc = 1:data.numberOfCarriers
         f[icc] = u[icc] - 0.0
        end

    else
        for icc in [iphin, iphip] 

            if data.recombinationOn == true   
                # radiative recombination
                kernelRadiative = data.recombinationRadiative[ireg]
                
                # Auger recombination
                kernelAuger     = (data.recombinationAuger[iphin, ireg] * n + data.recombinationAuger[iphip, ireg] * p)
                
                # SRH recombination
                kernelSRH       = 1.0 / (  data.recombinationSRHLifetime[iphip, ireg] * (n + data.recombinationSRHTrapDensity[iphin, ireg]) + data.recombinationSRHLifetime[iphin, ireg] * (p + data.recombinationSRHTrapDensity[iphip, ireg]) )
                
                # full recombination
                f[icc]          = q* data.chargeNumbers[icc]* (kernelRadiative + kernelAuger + kernelSRH)*  excessCarrierDensTerm  - q * data.chargeNumbers[icc] * generation(data, ireg)

            else
                f[icc]          = 0.0 #- q * data.chargeNumbers[icc] * generation(data, ireg) 
            end
        end

        for icc in iphip+1:data.numberOfCarriers
            f[icc]              = u[icc] - 0.0
        end
    
    end
    
    f[ipsi]                     = - q * data.λ1 * f[ipsi]
end


function reactionGraded!(f, u, node, data)

    # indices
    iphin                 = 1
    iphip                 = 2
    ipsi                  = data.numberOfCarriers + 1            # final index for electrostatic potential
    ireg                  = node.region
    inode                 = node.index

    n                     = computeDensities(u, data, inode, ireg, iphin, ipsi, true)  # true for interior region
    p                     = computeDensities(u, data, inode, ireg, iphip, ipsi, true) 
    exponentialTerm       = exp((q *u[iphin] - q  * u[iphip]) / (kB*data.temperature))
    excessCarrierDensTerm = n*p * (1.0 - exponentialTerm)

    # rhs of NLP (charge density)
    for icc = 1:data.numberOfCarriers
        eta     = etaFunction(u, node, data, icc, ipsi) 
        #eta     = data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + (data.bandEdgeEnergy[ireg] + data.bandEdgeEnergyNode[icc, inode]) / q )
    
        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * (data.doping[icc, node.region] + data.dopingNode[icc, node.index])  # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * (data.densityOfStates[icc, node.region] + data.densityOfStatesNode[icc, node.index]) * data.F[icc](eta)   # add charge carrier

    end

    # rhs of continuity equations for electron and holes (bipolar reaction)
    
    if data.inEquilibrium == true 

        for icc = 1:data.numberOfCarriers
         f[icc] = u[icc] - 0.0
        end

    else
        for icc in [iphin, iphip] 

            if data.recombinationOn == true   
                
                # radiative recombination
                kernelRadiative = (data.recombinationRadiative[ireg] + data.recombinationRadiativeNode[inode])
                
                # Auger recombination
                kernelAuger     = (data.recombinationAuger[iphin, ireg] * n + data.recombinationAuger[iphip, ireg] * p)
                
                # SRH recombination
                kernelSRH       = 1.0 / (  data.recombinationSRHLifetime[iphip, ireg] * (n + data.recombinationSRHTrapDensity[iphin, ireg]) + data.recombinationSRHLifetime[iphin, ireg] * (p + data.recombinationSRHTrapDensity[iphip, ireg]) )
                
                # full recombination
                f[icc]          = q* data.chargeNumbers[icc]* (kernelRadiative + kernelAuger + kernelSRH)*  excessCarrierDensTerm  - q * data.chargeNumbers[icc] * generation(data, ireg)

            else
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
$(SIGNATURES)

The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have ``f[c_i] =  z_i  q ∂_t c_i`` and for the electrostatic potential ``f[ψ] = 0``.
[NOT TESTED.]
"""

function storage!(f, u, node, data)
    ipsi       = data.numberOfCarriers + 1
    
    for icc = 1:data.numberOfCarriers

        eta    = etaFunction(u, node, data, icc, ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = data.chargeNumbers[icc] * data.densityOfStates[icc, node.region] * data.F[icc](eta)
    end

    f[ipsi] = 0.0
end

"""
$(SIGNATURES)

Compute trap densities. Currently, only done for the Boltzmann statistics.

"""
function trapDensity(icc, ireg, data, Et) 
    data.densityOfStates[icc, ireg] * exp( data.chargeNumbers[icc] * (data.bandEdgeEnergy[icc, ireg] - Et) / (kB * data.temperature)) # need to subtract Eref
end


"""
$(SIGNATURES)

The classical Scharfetter-Gummel flux scheme.

"""

function ScharfetterGummel!(f, u, edge, data)
    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    
    dpsi     = ul[ipsi] - uk[ipsi]
    f[ipsi]  = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)
        etal               = etaFunction(ul, edge, data, icc, ipsi) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        nodel              = edge.node[2]
        nodek              = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        bp, bm             = fbernoulli_pm(data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT)
        f[icc]             = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
    
    end


end

"""
$(SIGNATURES)

The classical Scharfetter-Gummel flux scheme for space-dependent DOS and band-edge energies.

"""

function ScharfetterGummelGraded!(f, u, edge, data)
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

        j0                 =  data.chargeNumbers[icc] * q * data.UT

        etak               = data.chargeNumbers[icc] / data.UT * ( (uk[icc] - uk[ipsi]) + (data.bandEdgeEnergyNode[icc, nodek] + data.bandEdgeEnergy[icc, ireg]) / q )
        etal               = data.chargeNumbers[icc] / data.UT * ( (ul[icc] - ul[ipsi]) + (data.bandEdgeEnergyNode[icc, nodel] + data.bandEdgeEnergy[icc, ireg]) / q )


        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]
        mobilityl          = (data.mobility[icc, ireg] + data.mobilityNode[icc, nodel])
        mobilityk          = (data.mobility[icc, ireg] + data.mobilityNode[icc, nodek])
        densityOfStatesl   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodel])
        densityOfStatesk   = (data.densityOfStates[icc, ireg] + data.densityOfStatesNode[icc,nodek])
        
        if data.densityOfStatesNode[icc, nodel] ≈ 0.0 || data.densityOfStatesNode[icc, nodek] ≈ 0.0
            bp, bm         = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT ) 
        else
            bp, bm         = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / data.UT - (log(data.densityOfStatesNode[icc, nodel]) -log(data.densityOfStatesNode[icc,nodek])) ) 
        end

        f[icc]             = - j0  * ( bm  * mobilityl * densityOfStatesl * data.F[icc](etal) - bp * mobilityk * densityOfStatesk * data.F[icc](etak) )

    end

end


"""
$(SIGNATURES)

The Sedan flux scheme.

"""

function Sedan!(f, u, edge, data)
    
    uk       = viewK(edge, u)
    ul       = viewL(edge, u)
    
    ipsi     = data.numberOfCarriers + 1
    ireg     = edge.region
    
    dpsi     = ul[ipsi] - uk[ipsi]
    f[ipsi]  = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi) # calls etaFunction(u, edge, data, icc, ipsi)
        etal               = etaFunction(ul, edge, data, icc, ipsi) # calls etaFunction(u, edge, data, icc, ipsi)

        nodel              = edge.node[2]
        nodek              = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        Q                  = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm             = fbernoulli_pm(Q)

        f[icc] = - j0 * ( bm * data.F[icc](etal) - bp * data.F[icc](etak) )
end

end

"""
$(SIGNATURES)

The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is 
used for the regularization of the removable singularity.

"""

function diffusionEnhanced!(f, u, edge, data)
    tolReg  = 1.0e-13;
    
    uk      = viewK(edge, u)
    ul      = viewL(edge, u)
    
    ipsi    = data.numberOfCarriers + 1
    ireg    = edge.region
    
    dpsi    = ul[ipsi] - uk[ipsi]
    f[ipsi] = - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers

        j0                 =  data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak               = etaFunction(uk, edge, data, icc, ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
        etal               = etaFunction(ul, edge, data, icc, ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)

        if abs( (etal - etak)/(etak + etal) ) > tolReg
            g  = (etal - etak ) / ( log(data.F[icc](etal)) - log(data.F[icc](etak)) )
        else
            # regularization idea coming from Pietra-Jüngel scheme
            gk = exp(etak) / data.F[icc](etak)
            gl = exp(etal) / data.F[icc](etal)
            g  = 0.5 * ( gk + gl )
        end

        nodel              = edge.node[2]
        nodek              = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        bp, bm             = fbernoulli_pm(data.chargeNumbers[icc] * (dpsi - bandEdgeDifference / q) / (data.UT * g))
        f[icc]             = - j0 * g * (  bm * data.F[icc](etal) - bp * data.F[icc](etak))
    end

end

"""
$(SIGNATURES)

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
    
    ipsi          = data.numberOfCarriers + 1
    ireg          = edge.region
    
    dpsi          = ul[ipsi] - uk[ipsi]
    f[ipsi]        =  - data.dielectricConstant[ireg] * ε0 * dpsi
    
    if data.inEquilibrium == true # return zero flux in equilibrium
        return
    end
    
    for icc = 1:data.numberOfCarriers 

        j0                  = data.chargeNumbers[icc] * q * data.mobility[icc, ireg] * data.UT * data.densityOfStates[icc, ireg]

        etak                = etaFunction(uk, edge, data, icc, ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
        etal                = etaFunction(ul, edge, data, icc, ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)

        # use Sedan flux as starting guess

        nodel               = edge.node[2]
        nodek               = edge.node[1]
        bandEdgeDifference  = data.bandEdgeEnergyNode[icc, nodel] - data.bandEdgeEnergyNode[icc, nodek]

        Q                   = data.chargeNumbers[icc]*( (dpsi - bandEdgeDifference/q) /data.UT) + (etal - etak) - log(data.F[icc](etal)) + log(data.F[icc](etak) )
        bp, bm              = fbernoulli_pm(Q)
        jInitial            = ( bm * data.F[icc](etal)  - bp * data.F[icc](etak))

        implicitEq(j::Real) = (fbernoulli_pm(data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q) /data.UT) - data.γ * j)[2] * exp(etal) - fbernoulli_pm(data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q) /data.UT) - data.γ * j)[1] * exp(etak)) - j

        delta               = 1.0e-18 + 1.0e-14 * abs(value(jInitial))

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

            it       = it + 1
            damp     = min(damp * 1.2, 1.0)
        end
        f[icc]       = - j0 * jInitial
    end
end


"""

$(SIGNATURES)

The argument of the distribution function for floats.

    z / UT  * ( (phi - psi) + E / q ).

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

$(SIGNATURES)

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

$(SIGNATURES)

For given potentials in vector form, compute corresponding vectorized densities.
(Caution: DID NOT TAKE JUMPING PARAMETERS INTO ACCOUNT!!!!)
[NOT TESTED FOR MULTIDIMENSIONS]
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

For given solution in vector form, compute corresponding vectorized band-edge energies and fermi level.
(Caution: DID NOT TAKE JUMPING PARAMETERS INTO ACCOUNT!!!!)
[NOT TESTED FOR MULTIDIMENSIONS]

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

$(SIGNATURES)

Compute the electro-neutral solution for the Boltzmann approximation. 
It is obtained by setting the left-hand side in
the Poisson equation equal to zero and solving for \\psi.
The charge carriers may obey different statitics functions.

"""
function electroNeutralSolution!(data, grid; Newton=false)

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

$(SIGNATURES)

Compute the charge density, i.e. the right-hand side of Poisson's equation.

"""
function chargeDensity(psi0, phi, UT, EVector, chargeNumbers, dopingVector, dosVector, FVector)
    # https://stackoverflow.com/questions/45667291/how-to-apply-one-argument-to-arrayfunction-1-element-wise-smartly-in-julia
    sum(-chargeNumbers .* dopingVector) + sum(chargeNumbers .* dosVector .* (etaFunction(psi0, phi, UT, EVector, chargeNumbers) .|> FVector))
end

"""

$(SIGNATURES)

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
