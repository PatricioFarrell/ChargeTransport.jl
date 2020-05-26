"""
$(TYPEDEF)

Struct holding physical data for drift-diffusion simulation of semiconductor device.
If there are ``N`` number of species, it is assumed that the first ``N-1``ones
correspond to the charge carriers and the final one to the electrostatic potential.

$(TYPEDFIELDS)
"""

mutable struct DDFermiData <: VoronoiFVM.AbstractData

    # integer numbers
    numberOfNodes               ::  Int64
    numberOfRegions             ::  Int64
    numberOfBoundaryRegions     ::  Int64
    numberOfSpecies             ::  Int64

    # distribution (Boltzmann, Blakemore, Fermi-Dirac etc.)
    F                           ::  Function

    # real numbers
    temperature                 ::  Real
    UT                          ::  Real

    # number of boundary regions
    contactVoltage              ::  Array{Real,1}

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
    dielectricConstant          ::  Array{Real,1}
#da: warum recombination Radiative nur ein dimensional? -> Rekombinationen anschauen.
    recombinationRadiative      ::  Array{Real,1}
    electronSpinRelaxationTime  ::  Array{Real,1}
    holeSpinRelaxationTime      ::  Array{Real,1}
    recombinationDirect         ::  Array{Real,1}
    generationEmittedLight      ::  Array{Real,1}
    generationPrefactor         ::  Array{Real,1}
    generationAbsorption        ::  Array{Real,1}

    # number of nodes x number of carriers
    dopingNode                  ::  Array{Real,2}
    densityOfStatesNode         ::  Array{Real,2}   # still needs to be implemented
    bandEdgeEnergyNode          ::  Array{Real,2}   # still needs to be implemented

    # standard constructor
    # DDFermiData(... all args ...) = new(... all args ...)

end

function emptyFunction()
end

"""

$(SIGNATURES)

Simplified constructors for DDFermiData which takes only the
number of regions, number of boundary regions and the number
of charge carriers as input.

"""

function DDFermiData(numberOfNodes::Int64, numberOfRegions=3::Int64, numberOfBoundaryRegions=2::Int64, numberOfSpecies=3::Int64)
    DDFermiData(

    # integer numbers
    numberOfNodes,
    numberOfRegions,
    numberOfBoundaryRegions,
    numberOfSpecies,

    # functions
    emptyFunction,                                                  # distribution

    # real numbers
    300 * K,                                                        # temperature
    (kB * 300 * K ) / q,                                            # thermal voltage

    # number of boundary regions
    Array{Real,1}(undef,numberOfBoundaryRegions),                   # contactVoltage

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
    Array{Real,1}(undef,numberOfRegions),                           # dielectricConstant
    Array{Real,1}(undef,numberOfRegions),                           # recombinationRadiative
    Array{Real,1}(undef,numberOfRegions),                           # electronSpinRelaxationTime
    Array{Real,1}(undef,numberOfRegions),                           # holeSpinRelaxationTime
    Array{Real,1}(undef,numberOfRegions),                           # recombinationDirect
    Array{Real,1}(undef,numberOfRegions),                           # generationEmittedLight
    Array{Real,1}(undef,numberOfRegions),                           # generationPrefactor
    Array{Real,1}(undef,numberOfRegions),                           # generationAbsorption

    # number of nodes x number of carriers
    spzeros(Float64,numberOfNodes,numberOfSpecies-1),               # dopingNode
    spzeros(Float64,numberOfNodes,numberOfSpecies-1),               # densityOfStatesNode
    spzeros(Float64,numberOfNodes,numberOfSpecies-1)                # bandEdgeEnergyNode
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
function etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
    E  = data.bandEdgeEnergy[node.region,icc] + data.bandEdgeEnergyNode[node.index,icc]
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for boundary nodes:
    z / UT  * ( (phi_at_boundary - psi) + E / q ).
"""

function etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
    E  = data.bBandEdgeEnergy[bnode.region,icc] + data.bandEdgeEnergyNode[bnode.index,icc]
    data.chargeNumbers[icc] / data.UT * ( (data.contactVoltage[bnode.region]- u[ipsi]) + E / q )
end

"""

$(SIGNATURES)

The argument of the distribution function for edges:

    z / UT  * ( (phi_at_edge - psi) + E / q ).

"""

function etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
    E  = data.bandEdgeEnergy[edge.region,icc] + data.bandEdgeEnergyNode[edge.icell,icc] #icell: Number of discretization cell the edge is invoked from
    data.chargeNumbers[icc] / data.UT * ( (u[icc] - u[ipsi]) + E / q )
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
    ipsi = data.numberOfSpecies             # final index for electrostatic potential

    for icc = 1:data.numberOfSpecies - 1

        eta = etaFunction(u,bnode,data,icc,ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)

        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.bDoping[bnode.region,icc]                            # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.bDensityOfStates[bnode.region,icc] * data.F(eta)     # add charge carrier

        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0

    end
    f[ipsi] = -1/α *  q * f[ipsi]

end
"""
$(SIGNATURES)

Boundary reaction term for densities.
"""

function breactionDensities!(f,u,bnode,data)

    # parameters
    α    = 1.0/VoronoiFVM.Dirichlet         # tiny penalty value
    ipsi = data.numberOfSpecies             # final index for electrostatic potential

    for icc = 1:data.numberOfSpecies - 1

        eta = etaFunction(u,bnode,data,icc,ipsi) # calls etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)

        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.bDoping[bnode.region,icc]                            # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.bDensityOfStates[bnode.region,icc] * data.F(eta)

        # boundary conditions for charge carriers are set in main program
        f[icc]  = 0.0# data.bDensityOfStates[bnode.region,icc] * data.F(eta) - data.chargeNumbers[icc] * u[icc]

    end
    f[ipsi] = -1/α *  q * f[ipsi]

end

"""
$(SIGNATURES)

Sets up the right-hand sides. Assuming a bipolar semiconductor
the right-hand side for the electrostatic potential becomes

  ``f[ψ]  = - q ((p - N_a) - (n - N_d) ) = - q  \\sum  z_i  (c_i - N_i) ``

and the right-hand sides for the charge carriers yields

``f[c_i] =  z_i  q  R ``

for a charge number ``z_i`` and all charge carriers ``c_i``.
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
    ipsi = data.numberOfSpecies             # final index for electrostatic potential

    for icc = 1:data.numberOfSpecies - 1

        eta = etaFunction(u,node,data,icc,ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)

        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.doping[node.region,icc]                          # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * data.densityOfStates[node.region,icc] * data.F(eta)   # add charge carrier

        ## add different recombination kernels r(n,p)
        for ireg = 1:data.numberOfRegions

            # radiative recombination
            f[icc] = data.recombinationRadiative[ireg]

            # Auger recombination
            f[icc] = f[icc] + sum(data.recombinationAuger[ireg,:] .* u[1:end-1])

            # SRH recombination
            # da: TrapDensity und LifeTime vertauscht.
            #vllt sowas wie data.recombinationSRHActive ?
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

"""
$(SIGNATURES)

Reaction term for densities. For simplicity, the recombination is neglected.
"""

function reactionDensities!(f,u,node,data)

    # parameters
    ipsi = data.numberOfSpecies             # final index for electrostatic potential

    for icc = 1:data.numberOfSpecies - 1

        f[ipsi] = f[ipsi] - data.chargeNumbers[icc] * data.doping[node.region,icc]                          # subtract doping
        f[ipsi] = f[ipsi] + data.chargeNumbers[icc] * u[icc]   # add charge carrier
        f[icc]  = 0.0

    end
    f[ipsi] = - q * f[ipsi]
end

"""
$(SIGNATURES)

The storage term for time-dependent problems.
Currently, for the time-dependent current densities the implicit Euler scheme is used.
Hence, we have ``f[c_i] =  z_i  q ∂_t c_i`` and for the electrostatic potential ``f[ψ] = 0``.
"""

function storage!(f, u, node, data)
    ipsi = data.numberOfSpecies

    for icc = 1:data.numberOfSpecies - 1
        eta = etaFunction(u,node,data,icc,ipsi) # calls etaFunction(u,node::VoronoiFVM.Node,data,icc,ipsi)
        f[icc] = data.chargeNumbers[icc] * data.densityOfStates[node.region, icc] * data.F(eta)
    end
    f[ipsi] = 0.0
end
"""
$(SIGNATURES)

The classical Scharfetter-Gummel flux scheme for densities.

"""

function ScharfetterGummelDensities!(f, u, edge, data)
    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    dpsi     = ul[ipsi]- uk[ipsi]

    f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

    for icc = 1:data.numberOfSpecies-1

        j0    =  - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * data.UT

        nodel = edge.node[2]
        nodek = edge.node[1]

        bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

        bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference/q)/ data.UT)
        f[icc] = - data.chargeNumbers[icc] * j0 * ( bm * ul[icc] - bp * uk[icc] )

    end

end


"""
$(SIGNATURES)

The classical Scharfetter-Gummel flux scheme.

"""

function ScharfetterGummel!(f, u, edge, data)
    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    dpsi     = ul[ipsi]- uk[ipsi]

    f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

    for icc = 1:data.numberOfSpecies-1

        j0    =  - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * data.UT * data.densityOfStates[ireg,icc]

        etak  = etaFunction(uk,edge,data,icc,ipsi) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)
        etal  = etaFunction(ul,edge,data,icc,ipsi) # calls etaFunction(u, edge::VoronoiFVM.Edge, data, icc, ipsi)

        nodel = edge.node[2]
        nodek = edge.node[1]

        bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

        bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference/q)/ data.UT)
        f[icc] = - data.chargeNumbers[icc] * j0 * ( bm * data.F(etal) - bp * data.F(etak) )

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

    dpsi     = ul[ipsi]- uk[ipsi]
    f[ipsi]  = - data.dielectricConstant[ireg] * ε0 * dpsi

    for icc = 1:data.numberOfSpecies-1

        j0       = - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * data.UT * data.densityOfStates[ireg,icc]

        etak = etaFunction(uk,edge,data,icc,ipsi) # calls etaFunction(u, edge, data, icc, ipsi)
        etal = etaFunction(ul,edge,data,icc,ipsi) # calls etaFunction(u, edge, data, icc, ipsi)

        nodel = edge.node[2]; nodek = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

        Q = data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q)/ data.UT) + (etal - etak) - log(data.F(etal)) + log(data.F(etak))
        bp, bm = fbernoulli_pm(Q)

        f[icc] = data.chargeNumbers[icc] * j0 * ( bm * data.F(etal) - bp * data.F(etak) )
    end

end

"""
$(SIGNATURES)

The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is used for regularization of removable singularity.

"""

function diffusionEnhanced!(f, u, edge, data)
    tolRegularisation = 1.0e-13;

    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    dpsi = ul[ipsi]- uk[ipsi]
    f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

    for icc = 1:data.numberOfSpecies-1

        j0   = - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * data.UT * data.densityOfStates[ireg,icc]

        etak = etaFunction(uk,edge,data,icc,ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
        etal = etaFunction(ul,edge,data,icc,ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)

        if abs( (etal-etak)/ (etak+etal)) > tolRegularisation
            g = (etal - etak ) / ( log(data.F(etal)) - log(data.F(etak)) )
        else # regularization idea coming from Pietra-Jüngel scheme
            gk = exp(etak)/data.F(etak)
            gl = exp(etal)/data.F(etal)
            g  = 0.5 * ( gk + gl )
        end

        nodel = edge.node[2]; nodek = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

        bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference/q)/ (data.UT * g) )
        f[icc] = - data.chargeNumbers[icc] * j0 * g * (  bm * data.F(etal) - bp * data.F(etak))
    end

end

"""
$(SIGNATURES)

The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation which arises when considering the generalized Scharfetter-Gummel scheme in case of Blakemore statistics.
Hence, it should be exclusively worked with, when considering the Blakemore distribution.

"""

function KopruckiGaertner!(f, u, edge, data)
    gamma = 0.27        # from Blakemore distribution
    max_iteration = 200 # for Newton solver
    it = 0              # number of iterations (newton)
    damp = 0.1          # damping factor

    uk  = viewK(edge, u)
    ul  = viewL(edge, u)

    ipsi = data.numberOfSpecies
    ireg = edge.region

    dpsi     = ul[ipsi]- uk[ipsi]
    f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

    for icc = 1:data.numberOfSpecies-1

        j0   = - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * data.UT * data.densityOfStates[ireg,icc]

        etak = etaFunction(uk,edge,data,icc,ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
        etal = etaFunction(ul,edge,data,icc,ipsi) # calls etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)

        # use Sedan flux as starting guess
        Q = data.chargeNumbers[icc] * dpsi/ data.UT + (etal - etak) - log(data.F(etal)) + log(data.F(etak))
        bp, bm = fbernoulli_pm(Q)
        jInitial =  ( bm * data.F(etal)  - bp * data.F(etak))

        implicitEquation(j::Real) =  (fbernoulli_pm(data.chargeNumbers[icc] * (dpsi / data.UT) - gamma*j )[2] * exp(etal) - fbernoulli_pm(data.chargeNumbers[icc] * (dpsi/ data.UT) - gamma*j )[1] * exp(etak)) - j

        delta = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
        while (it < max_iteration)
            Fval  = implicitEquation(jInitial)
            dFval = ForwardDiff.derivative(implicitEquation,jInitial)
            if isnan( value(dFval) ) || value( abs(dFval) ) < delta
                @show value(jInitial), value(Fval), value(dFval)
                error("singular derivative")
            end
            update = Fval/dFval
            jInitial = jInitial - damp * update
            if abs(update) < delta
                break
            end
            it = it + 1
            damp = min(damp*1.2, 1.0)
        end
        f[icc] =   data.chargeNumbers[icc] * j0 * jInitial
    end
end

"""

$(SIGNATURES)

For given potentials, compute corresponding densities.

"""
function calculateDensities(grid, sys, U0)
    dddata        = data(sys)
    ipsi          = dddata.numberOfSpecies

    coord         = grid[Coordinates]
    bfaceregions  = grid[BFaceRegions]
    bfacenodes    = grid[BFaceNodes]
    cellregions   = grid[CellRegions]
    cellnodes     = grid[CellNodes]
    numberOfCoord = length(coord)

    # densities extrahieren in eigene Funktion.
    densities = Array{Real,2}(undef, dddata.numberOfSpecies-1, length(coord))

    for icc = 1:dddata.numberOfSpecies-1
        # Möglichkeit densities in einer zeile zu implementieren? (Über etaF()-> multiple dispatch)
        # --> versuche Informationen über coord zu bekommen!

        E = dddata.bBandEdgeEnergy[bfaceregions[1],icc] + dddata.bandEdgeEnergyNode[bfacenodes[1],icc]
        eta = dddata.chargeNumbers[icc] / dddata.UT * ( (dddata.contactVoltage[1]- U0[ipsi, 1]) + E / q )
        densities[icc, 1] = dddata.bDensityOfStates[1, icc] * dddata.F(eta)

        for i = 2:numberOfCoord-1
            # frage: wie am besten etaF() benutzen? -> etaF() hat als Eingeparameter etwas vom Typ VoronoiFVM.Node
            E   = dddata.bandEdgeEnergy[cellregions[i], icc] + dddata.bandEdgeEnergyNode[i, icc]
            eta = dddata.chargeNumbers[icc] / dddata.UT * ( (U0[icc, i]- U0[ipsi, i]) + E / q )
            densities[icc, i] = dddata.densityOfStates[cellregions[i], icc] * dddata.F(eta)

        end
        E = dddata.bBandEdgeEnergy[bfaceregions[2],icc] + dddata.bandEdgeEnergyNode[bfacenodes[2],icc]
        eta = dddata.chargeNumbers[icc] / dddata.UT * ( (dddata.contactVoltage[2]- U0[ipsi, numberOfCoord]) + E / q )
        densities[icc, numberOfCoord] = dddata.bDensityOfStates[2, icc] * dddata.F(eta)

    end
    return densities
end

"""

$(SIGNATURES)

Compute the electro-neutral solution, assuming the Boltzmann approximation. It is obtained by setting the left-hand side in
the Poisson equation equal to zero and solving for \\psi.

"""

function electroNeutralSolutionBoltzmann(grid::ExtendableGrid,data::DDFermiData)

    if data.numberOfSpecies-1 != 2
        error("The electroneutral solution is only implemented for two species!")
    end

    # region independent parameters
    iphin = 1;
    iphip = 2;

    # initialize zero vector
    coord        = grid[Coordinates]
    bfacenodes   = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    psi0         = zeros(length(coord))

    # boundary values
    for i = 1:length(bfacenodes)
        # boundary index
        ibreg = bfaceregions[i]

        # boundary region specific data
        Ec = data.bBandEdgeEnergy[ibreg,iphin]
        Ev = data.bBandEdgeEnergy[ibreg,iphip]
        Nc = data.bDensityOfStates[ibreg,iphin]
        Nv = data.bDensityOfStates[ibreg,iphip]
        Ni = sqrt( Nc*Nv*exp(-(Ec-Ev)/(kB*data.temperature)) )
        C  = -data.chargeNumbers[iphin] * data.bDoping[ibreg,iphin] -data.chargeNumbers[iphip] * data.bDoping[ibreg,iphip]

        # set boundary values for electroneutral potential
        psi0[bfacenodes[i]] = (Ec+Ev)/(2q) - 0.5*data.UT*log(Nc/Nv) + data.UT*asinh(C/(2*Ni))
    end

    # interior values
    cellregions = grid[CellRegions]

    for i=1:length(cellregions)-1

        # interior index
        ireg      = cellregions[i]
        ireg_next = cellregions[i+1]

        # interior region specific data
        Ec = (data.bandEdgeEnergy[ireg,iphin]  + data.bandEdgeEnergy[ireg_next,iphin] ) / 2 + data.bandEdgeEnergyNode[i,iphin]
        Ev = (data.bandEdgeEnergy[ireg,iphip]  + data.bandEdgeEnergy[ireg_next,iphip] ) / 2 + data.bandEdgeEnergyNode[i,iphip]
        Nc = (data.densityOfStates[ireg,iphin] + data.densityOfStates[ireg_next,iphin]) / 2
        Nv = (data.densityOfStates[ireg,iphip] + data.densityOfStates[ireg_next,iphip]) / 2
        Ni = sqrt( Nc*Nv*exp(-(Ec-Ev)/(kB*data.temperature)) )
        C  = -data.chargeNumbers[iphin] * (data.doping[ireg,iphin]+data.doping[ireg_next,iphin])/2 -
              data.chargeNumbers[iphip] * (data.doping[ireg,iphip]+data.doping[ireg_next,iphip])/2

        # set interior values for electroneutral potential
        psi0[i+1] = (Ec+Ev)/(2q) - 0.5*data.UT*log(Nc/Nv) + data.UT*asinh(C/(2*Ni))
    end

    psi0
end

"""

$(SIGNATURES)

Find the equilibrium solution for the electrostatic potential

"""

function solveEquilibriumBoltzmann!(solution, initialGuess, data, grid, control, dense)

    if !(data.F == DDFermi.Boltzmann) # if F != Boltzmann, find equilibrium solution for Boltzmann
        num_cellregions = grid[NumCellRegions]
        num_bfaceregions = grid[NumBFaceRegions]
        species  = 1:data.numberOfSpecies
        regions  = 1:num_cellregions
        bregions = 1:num_bfaceregions

        # save and set new values (careful with aliasing of arrays!)
        saveDistribution    = data.F
        saveContactVoltage  = copy(data.contactVoltage)             # copy() avoids aliasing
        data.F              = Boltzmann
        data.contactVoltage = zeros(size(data.contactVoltage)) * V

        # initializing physics environment with the Boltzmann approximation as distribution function
        physicsBoltzmann = VoronoiFVM.Physics(
            data        = data,
            num_species = data.numberOfSpecies,
            flux        = DDFermi.ScharfetterGummel!,
            reaction    = DDFermi.reaction!,
            breaction   = DDFermi.breaction!
        )

        if dense
            sysBoltzmann = VoronoiFVM.System(grid, physicsBoltzmann, unknown_storage = :dense)
        else
            sysBoltzmann = VoronoiFVM.System(grid, physicsBoltzmann, unknown_storage = :sparse)
        end
        # enable all species in all regions
        for ispecies in species
            enable_species!(sysBoltzmann, ispecies, regions)
        end
        for icc in species[1:end-1]
            for bregion in bregions
                sysBoltzmann.boundary_values[icc,  bregion] = data.contactVoltage[bregion]
                sysBoltzmann.boundary_factors[icc, bregion] = VoronoiFVM.Dirichlet
            end
        end
        solve!(solution, initialGuess, sysBoltzmann, control = control, tstep = Inf)
        initialGuess .= solution

        # switch back to the original data values
        data.F              = saveDistribution
        data.contactVoltage = saveContactVoltage

    else # if F = Boltzmann, don't do anything

        println("*** We compute with Boltzmann statistics anyway. ")

    end

end
