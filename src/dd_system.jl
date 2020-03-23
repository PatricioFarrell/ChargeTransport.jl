"""
$(TYPEDEF)

Struct holding physical data for drift-diffusion simulation of semiconductor device.
If there are ``N`` number of species, it is assumed that the first ``N-1``ones 
correspond to the charge carriers and the final one to the electrostatic potential. 

$(TYPEDFIELDS)
"""
mutable struct DDFermiData <: VoronoiFVM.AbstractData
    
    # integer numbers
    numberOfRegions             ::  Int64
    numberOfBoundaryRegions     ::  Int64
    numberOfSpecies             ::  Int64

    # distribution (Boltzmann, Blakemore, Fermi-Dirac etc.)
    F                           ::  Function                  

    # real numbers
    temperature                 ::  Real
    
    # number of boundary regions 
    contactVoltage              ::  Array{Real,1}
    bDopingClassical            ::  Array{Real,1}

    # number of carriers
    chargeNumbers               ::  Array{Real,1}

    # number of boundary regions x number of carriers
    bBandEdgeEnergy             ::  Array{Real,2}
    bDensityOfStates            ::  Array{Real,2}
    bDoping                     ::  Array{Real,2}
    
    # number of regions x number of carriers
    doping                      ::  Array{Real,2}
    densityOfStates             ::  Array{Real,2}
    bandEdgeEnergy              ::  Array{Real,2}
    mobility                    ::  Array{Real,2}
    recombinationSRHLifetime    ::  Array{Real,2}
    recombinationSRHTrapDensity ::  Array{Real,2}
    recombinationAuger          ::  Array{Real,2}

    # number of regions 
    # intrinsicDoping             ::  Array{Real,1}
    dopingClassical             ::  Array{Real,1}
    dielectricConstant          ::  Array{Real,1}
    recombinationRadiative      ::  Array{Real,1}
    electronSpinRelaxationTime  ::  Array{Real,1}
    holeSpinRelaxationTime      ::  Array{Real,1}
    recombinationDirect         ::  Array{Real,1}
    generationEmittedLight      ::  Array{Real,1}
    generationPrefactor         ::  Array{Real,1}
    generationAbsorption        ::  Array{Real,1}

    # standard constructor
    # DDFermiData(... all args ...) = new(... all args ...)

end

function nofunc()
end

"""

$(SIGNATURES)

Simplified constructors for DDFermiData which takes only the 
number of regions, number of boundary regions and the number 
of charge carriers as input.

"""
function DDFermiData(numberOfRegions=3::Int64, numberOfBoundaryRegions=2::Int64, numberOfSpecies=3::Int64)
    DDFermiData(

        # integer numbers
        numberOfRegions,
        numberOfBoundaryRegions,
        numberOfSpecies,

        # functions
        nofunc,                                                         # distribution

        # real numbers
        300 * K,                                                        # temperature 

        # number of boundary regions
        Array{Real,1}(undef,numberOfBoundaryRegions),                   # contactVoltage
        Array{Real,1}(undef,numberOfBoundaryRegions),                   # bDopingClassical

        # number of charge carriers = number of species - 1
        Array{Real,1}(undef,numberOfSpecies-1),                         # chargeNumbers

        # number of boundary regions x number of carriers
        Array{Real,2}(undef,numberOfBoundaryRegions,numberOfSpecies-1), # bBandEdgeEnergy
        Array{Real,2}(undef,numberOfBoundaryRegions,numberOfSpecies-1), # bDensityOfStates
        zeros(Float64,      numberOfBoundaryRegions,numberOfSpecies-1), # bDoping

        # number of regions x number of charge carriers 
        zeros(Float64,      numberOfRegions,numberOfSpecies-1),         # doping 
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # densityOfStates
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # bandEdgeEnergy 
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # mobility     
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # recombinationSRHLifetime
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # recombinationSRHTrapDensity
        Array{Real,2}(undef,numberOfRegions,numberOfSpecies-1),         # recombinationAuger      

        # number of regions
        # zeros(Float64,      numberOfRegions),                           # intrinsicDoping
        Array{Real,1}(undef,numberOfRegions),                           # dopingClassical 
        Array{Real,1}(undef,numberOfRegions),                           # dielectricConstant      
        Array{Real,1}(undef,numberOfRegions),                           # recombinationRadiative  
        Array{Real,1}(undef,numberOfRegions),                           # electronSpinRelaxationTime
        Array{Real,1}(undef,numberOfRegions),                           # holeSpinRelaxationTime  
        Array{Real,1}(undef,numberOfRegions),                           # recombinationDirect     
        Array{Real,1}(undef,numberOfRegions),                           # generationEmittedLight  
        Array{Real,1}(undef,numberOfRegions),                           # generationPrefactor     
        Array{Real,1}(undef,numberOfRegions),                           # generationAbsorption  
        
    )

end

function Base.show(io::IO, this::DDFermiData)
    for name in fieldnames(typeof(this))
        @printf("%30s = ",name)
        println(io,getfield(this,name))
    end
end

"""

$(SIGNATURES)

The argument of the distribution function for interior nodes:

    z / UT  * ( (phi - psi) + E / q ).

"""
function etaFunction(u,node::VoronoiFVM.Node{Float64,Int32},data,icc,ipsi)
    UT = (kB * data.temperature ) / q 
    data.chargeNumbers[icc] / UT * ( (u[icc] - u[ipsi]) + data.bandEdgeEnergy[node.region,icc] / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for boundary nodes:

    z / UT  * ( (phi_at_boundary - psi) + E / q ).

"""
function etaFunction(u,bnode::VoronoiFVM.BNode{Float64,Int32},data,icc,ipsi)
    UT = (kB * data.temperature ) / q 
    data.chargeNumbers[icc] / UT * ( (data.contactVoltage[bnode.region]- u[ipsi]) + data.bandEdgeEnergy[bnode.region,icc] / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for edges:

    z / UT  * ( (phi_at_boundary - psi) + E / q ).

"""
function etaFunction(u,edge::VoronoiFVM.Edge{Float64,Int32},data,icc,ipsi)
    UT = (kB * data.temperature ) / q 
    data.chargeNumbers[icc] / UT * ( (u[icc] - u[ipsi]) + data.bandEdgeEnergy[edge.region,icc] / q )
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
function breaction!(f,u,bnode,data)

    # parameters
    α    = 1.0/VoronoiFVM.Dirichlet         # tiny penalty value
    UT   = (kB * data.temperature ) / q     # thermal voltage
    ipsi = data.numberOfSpecies             # final index for electrostatic potential 

    for icc = 1:data.numberOfSpecies - 1 

        # eta = data.chargeNumbers[icc] / UT  * ( (data.contactVoltage[bnode.region] - u[ipsi]) + data.bandEdgeEnergy[bnode.region,icc] / q )
        eta = etaFunction(u,bnode,data,icc,ipsi)
        
        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.bDoping[bnode.region,icc]                            # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.bDensityOfStates[bnode.region,icc] * data.F(eta)     # add charge carrier

        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0   

    end

    f[ipsi] = -1/α *  q * f[ipsi]

end 




"""
$(SIGNATURES)

Sets up the right-hand sides. Assuming a bipolar semiconductor
the right-hand side for the electrostatic potential becomes

    f[ipsi]  = - q * ((p - N_a) - (n - N_d) ) = - q * sum { z_i * (c_i - N_i) }

and the right-hand sides for the charge carriers yields

    f[icc] =  z_i * q * R

for a charge number `z_i` and all charge carriers `icc`. 
The recombination includes radiative, Auger and Shockley-Read-Hall 
recombination.

Also semiconductor devices with more than two species are permitted.
However, the implementation of the Shockley-Read-Hall kernel might not easily 
generalize to arbitrary numbers of species.

Currently, it is done as follows:

    1 / ( sum(data.recombinationSRHTrapDensity[ireg,end:-1:1] .* (u[1:end-1] .+ data.recombinationSRHLifetime[ireg,1:end] ) ) )

It needs to be used carefully.

"""
function reaction!(f,u,node,data)

    # parameters
    UT   = (kB * data.temperature ) / q     # thermal voltage
    ipsi = data.numberOfSpecies             # final index for electrostatic potential

    # set intrinsic doping outside of charge carrier loop
    # f[ipsi] = data.intrinsicDoping[node.region]

    for icc = 1:data.numberOfSpecies - 1 
        
        eta = etaFunction(u,node,data,icc,ipsi)
        # eta = data.chargeNumbers[icc] / UT * ( (u[icc] - u[ipsi]) + data.bandEdgeEnergy[node.region,icc]/q )

        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.doping[node.region,icc]                          # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.densityOfStates[node.region,icc] * data.F(eta)   # add charge carrier

        ## add different recombination kernels r(n,p)
        for ireg = 1:data.numberOfRegions

            # radiative recombination
            f[icc] = data.recombinationRadiative[ireg]           
            
            # Auger recombination
            f[icc] = f[icc] + sum(data.recombinationAuger[ireg,:] .* u[1:end-1])        

            # SRH recombination
            f[icc] = f[icc] + 1 / ( sum(data.recombinationSRHTrapDensity[ireg,end:-1:1] .* (u[1:end-1] .+ data.recombinationSRHLifetime[ireg,1:end] ) ) ) 
        end
        
        # full recombination
        # note: typeof(vec .* vec) is Array so we compute (vec .* vec)[1]
        f[icc]  = + q * data.chargeNumbers[icc] * f[icc] * prod(u[1:end-1]) * ( 1 - prod( exp( (-data.chargeNumbers .* u[1:end-1])[1] ) ) )

        # try 
        #     println(f[icc].value)
        # catch
        #     println(f[icc])
        # end
        
    end

    f[ipsi] = - q * f[ipsi]

end


# """
# $(SIGNATURES)

# Like `breaction!` but with classical regionwise doping.

# """
# function breaction_classicalRegionwiseDoping!(f,u,bnode,data)

#     # tiny penalty value
#     α = 1.0/VoronoiFVM.Dirichlet     

#     # doping and values for psi at Dirichlet boundary interfaces
#     bDopingVector   = data.bDopingClassical   # [p-doped, n-doped]
#     contactVoltages = data.contactVoltage     # [p-doped, n-doped]

#     # final index for electrostatic potential
#     ipsi = data.numberOfSpecies
    
#     # set up boundary conditions via penalty method
#     f[ipsi] = doping(bnode.region, bDopingVector)

#     for icc = 1:data.numberOfSpecies - 1 
#         eta = data.chargeNumbers[icc] * (q * (contactVoltages[bnode.region] - u[ipsi]) + data.bandEdgeEnergy[bnode.region,icc] )/ (kB * data.temperature)

#         f[icc]  = 0.0
#         f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.bDensityOfStates[bnode.region,icc] * data.F(eta)
#     end

#     f[ipsi] = -1/α * q * f[ipsi]

# end 


# """
# $(SIGNATURES)

# Like `reaction!` but with classical regionwise doping.

# """
# function reaction_classicalRegionwiseDoping!(f,u,node,data)

#     # final index for electrostatic potential
#     ipsi = data.numberOfSpecies

#     # extract doping from data
#     dopingVector   = data.dopingClassical

#     # set up right-hand sides
#     f[ipsi] = doping(node.region, dopingVector)


#     for icc = 1:data.numberOfSpecies - 1 
        
#         eta = data.chargeNumbers[icc] * (q * (u[icc] - u[ipsi]) + data.bandEdgeEnergy[node.region,icc] )/ (kB * data.temperature)

#         f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.densityOfStates[node.region,icc] * data.F(eta)

#         for ireg = 1:data.numberOfRegions

#             ## add different recombination kernels r(n,p)

#             # radiative recombination
#             f[icc] = data.recombinationRadiative[ireg]           
            
#             # Auger recombination
#             f[icc] = f[icc] + sum(data.recombinationAuger[ireg,:] .* u[1:end-1])        

#             # SRH recombination
#             f[icc] = f[icc] + 1 / ( sum(data.recombinationSRHTrapDensity[ireg,end:-1:1] .* (u[1:end-1] .+ data.recombinationSRHLifetime[ireg,1:end] ) ) )

#         end
        
#         # full recombination
#         # note: typeof(vec .* vec) is Array so we compute (vec .* vec)[1]
#         f[icc]  = + q * data.chargeNumbers[icc] * f[icc] * prod(u[1:end-1]) * ( 1 - prod( exp( (-data.chargeNumbers .* u[1:end-1])[1] ) ) )
        
#     end

#     f[ipsi] = - q * f[ipsi]

# end

# function doping(ireg,dopingVector)
#     dopingVector[ireg]
# end


###########################################################################################################################
########                                       DIFFERENT FLUX DISCRETIZATIONS                                      ########
###########################################################################################################################
"""
$(SIGNATURES)

The classical Scharfetter-Gummel flux scheme.

"""
function ScharfetterGummel!(f, u, edge, data)
    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    UT   = (kB * data.temperature ) / q
    dpsi = ul[ipsi]- uk[ipsi]

    for icc = 1:data.numberOfSpecies-1
        
        j0   = data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * UT * data.densityOfStates[ireg,icc]

        f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

        # etak = data.chargeNumbers[icc] / UT * ( uk[icc]-uk[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
        # etal = data.chargeNumbers[icc] / UT * ( ul[icc]-ul[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
        etak = etaFunction(uk,edge,data,icc,ipsi) 
        etal = etaFunction(ul,edge,data,icc,ipsi) 

        bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * dpsi / UT)

        f[icc] = data.chargeNumbers[icc] * j0 * ( bp * data.F(etak) - bm * data.F(etal) )

        # general implementation of the two equations:
        #       f[iphin] = - j0N * ( bp * data.F(etaNl) - bm * data.F(etaNk) )
        #       f[iphip] =   j0P * ( bp * data.F(etaPk) - bm * data.F(etaPl) )

    end

end

"""
$(SIGNATURES)

The Sedan flux scheme.

"""
function Sedan!(f, u, edge, data)
    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    UT   = (kB * data.temperature ) / q
    dpsi = ul[ipsi]- uk[ipsi]

    for icc = 1:data.numberOfSpecies-1

        j0   = data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * UT * data.densityOfStates[ireg,icc]

        f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

        # etak = data.chargeNumbers[icc] / UT * ( uk[icc]-uk[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
        # etal = data.chargeNumbers[icc] / UT * ( ul[icc]-ul[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
        etak = etaFunction(uk,edge,data,icc,ipsi) 
        etal = etaFunction(ul,edge,data,icc,ipsi) 

        Q = data.chargeNumbers[icc] / UT * dpsi + (etal - etak) - log( data.F(etal) ) + log( data.F(etak))

        bp, bm = fbernoulli_pm( Q )

        f[icc] = data.chargeNumbers[icc] * j0 * ( bp * data.F(etak) - bm * data.F(etal) )

    end

end

"""
Compute the electro-neutral solution, assuming the Boltzmann approximation. It is obtained by setting the left-hand side in 
the Poisson equation equal to zero and solving for \\psi.
"""
function electroNeutralSolutionBoltzmann(grid::VoronoiFVM.AbstractGrid,data::DDFermiData)

    if data.numberOfSpecies-1 != 2
        error("The electroneutral solution is only implemented for two species!")
    end

    # region independent parameters
    iphin = 1; 
    iphip = 2;
    UT    = (kB * data.temperature ) / q
    
    # initialize zero vector
    psi0 = zeros(length(grid.coord))
    
    # boundary values
    for i=1:length(grid.bfacenodes)

        # boundary index
        ibreg = grid.bfaceregions[i]

        # boundary region specific data
        Ec = data.bBandEdgeEnergy[ibreg,iphin]   
        Ev = data.bBandEdgeEnergy[ibreg,iphip] 
        Nc = data.bDensityOfStates[ibreg,iphin]  
        Nv = data.bDensityOfStates[ibreg,iphip]
        Ni = sqrt( Nc*Nv*exp(-(Ec-Ev)/(kB*data.temperature)) )
        C  = -data.chargeNumbers[iphin] * data.bDoping[ibreg,iphin] -data.chargeNumbers[iphip] * data.bDoping[ibreg,iphip]

        # set boundary values for electroneutral potential
        psi0[grid.bfacenodes[i]] = (Ec+Ev)/(2q) - 0.5*UT*log(Nc/Nv) + UT*asinh(C/(2*Ni))

    end

    # interior values
    for i=1:length(grid.cellregions)-1

        # interior index
        ireg      = grid.cellregions[i]
        ireg_next = grid.cellregions[i+1]

        # interior region specific data
        Ec = (data.bandEdgeEnergy[ireg,iphin]  + data.bandEdgeEnergy[ireg_next,iphin] ) / 2
        Ev = (data.bandEdgeEnergy[ireg,iphip]  + data.bandEdgeEnergy[ireg_next,iphip] ) / 2 
        Nc = (data.densityOfStates[ireg,iphin] + data.densityOfStates[ireg_next,iphin]) / 2
        Nv = (data.densityOfStates[ireg,iphip] + data.densityOfStates[ireg_next,iphip]) / 2
        Ni = sqrt( Nc*Nv*exp(-(Ec-Ev)/(kB*data.temperature)) )
        C  = -data.chargeNumbers[iphin] * (data.doping[ireg,iphin]+data.doping[ireg_next,iphin])/2 -
              data.chargeNumbers[iphip] * (data.doping[ireg,iphip]+data.doping[ireg_next,iphip])/2

        # set interior values for electroneutral potential
        psi0[i+1] = (Ec+Ev)/(2q) - 0.5*UT*log(Nc/Nv) + UT*asinh(C/(2*Ni))

    end

    psi0
end