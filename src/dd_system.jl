    """
    $(TYPEDEF)

    Struct holding physical data for drift-diffusion simulation of semiconductor device.
    If there are ``N`` number of species, it is assumed that the first ``N-1``ones
    correspond to the charge carriers and the final one to the electrostatic potential.

    $(TYPEDFIELDS)
    """
    mutable struct DDFermiData <: VoronoiFVM.AbstractData
    # todo_da: U_T hier ebenfalls data abhängig machen und in main file definieren, um nicht ständig auszurechen?
    # pf: okay!

        # integer numbers
        numberOfNodes               ::  Int64
        numberOfRegions             ::  Int64
        numberOfBoundaryRegions     ::  Int64
        numberOfSpecies             ::  Int64

        # distribution (Boltzmann, Blakemore, Fermi-Dirac etc.)
        F                           ::  Function

        # real numbers
        temperature                 ::  Real

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
        UT = (kB * data.temperature ) / q
        E  = data.bandEdgeEnergy[node.region,icc] + data.bandEdgeEnergyNode[node.index,icc]
        data.chargeNumbers[icc] / UT * ( (u[icc] - u[ipsi]) + E / q )
    end

    """

    $(SIGNATURES)

    The argument of the distribution function for boundary nodes:

        z / UT  * ( (phi_at_boundary - psi) + E / q ).

    """
    function etaFunction(u,bnode::VoronoiFVM.BNode,data,icc,ipsi)
        UT = (kB * data.temperature ) / q
        E  = data.bandEdgeEnergy[bnode.region,icc] + data.bandEdgeEnergyNode[bnode.index,icc]
        data.chargeNumbers[icc] / UT * ( (data.contactVoltage[bnode.region]- u[ipsi]) + E / q )
    #todo_da: muss hier nicht bBandEdgeEnergy[bnode.region, icc]?
    end

    """

    $(SIGNATURES)

    The argument of the distribution function for edges:

        z / UT  * ( (phi_at_edge - psi) + E / q ).

    """

    function etaFunction(u,edge::VoronoiFVM.Edge,data,icc,ipsi)
        UT = (kB * data.temperature ) / q
        E  = data.bandEdgeEnergy[edge.region,icc] + data.bandEdgeEnergyNode[edge.icell,icc]
    # todo_da: Kommentar: icell: Number of discretization cell the edge is invoked from
        data.chargeNumbers[icc] / UT * ( (u[icc] - u[ipsi]) + E / q )
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


    """
    $(SIGNATURES)

    The classical Scharfetter-Gummel flux scheme.

    """
    # todo_da:
    # -  j0 mit (-1) multiplizieren & f[icc] ebenfalls (für Konsistenz mit .tex Dateien)
    # -  f[ipsi] aus der Schleife nehmen (icc hat keinen Einfluss auf f[ipsi])
    # - anfangsbuchstabe klein, weil funktion        -

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

            nodel = edge.node[2]
            nodek = edge.node[1]

    # FRAGE: edge.node[2] = xL und edge.node[1] = xK?
    # todo_da: Ja, ergibt Sinn oder? Denn viewK(): "Solution view on first edge node" und mit neuer Schreibweise wird das auch deutlicher.
            bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

    # todo_da: bandEdgeDifference muss durch q dividiert werden
            bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference/q)/ UT)

    # todo_da: Kommentar: weiteres multiplizieren mit chargeNumber verursacht durch Kontinuitätsgleichung
            f[icc] = data.chargeNumbers[icc] * j0 * ( bp * data.F(etak) - bm * data.F(etal) )

            # general implementation of the two equations:
            #       f[iphin] = - j0N * ( bp * data.F(etaNl) - bm * data.F(etaNk) )
            #       f[iphip] =   j0P * ( bp * data.F(etaPk) - bm * data.F(etaPl) )
    # todo_da: oben genannte darstellung galt in meiner Implementation für:
    #    j0N = - data.q * data.muN * data.UT * data.Nc
    #    j0P =  data.q * data.muP * data.UT * data.Nv

        end

    end

    """
    $(SIGNATURES)

    The Sedan flux scheme.

    """
    # todo_da:
    # -  j0 mit (-1) multiplizieren & f[icc] ebenfalls (für Konsistenz mit .tex Dateien)
    # -  f[ipsi] aus der Schleife nehmen (icc hat keinen Einfluss auf f[ipsi])
    # - anfangsbuchstabe klein, weil funktion
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

            #todo_da:
            nodel = edge.node[2]; nodek = edge.node[1]
            bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]
            Q = data.chargeNumbers[icc] * ((dpsi - bandEdgeDifference/q)/ UT) + (etal - etak) - log(data.F(etal)) + log(data.F(etak))

            # Patricio:
            #Q = data.chargeNumbers[icc] / UT * dpsi + (etal - etak) - log( data.F(etal) ) + log( data.F(etak))
            bp, bm = fbernoulli_pm(Q)

            f[icc] = data.chargeNumbers[icc] * j0 * ( bp * data.F(etak) - bm * data.F(etal) )
        end

    end

    """
    $(SIGNATURES)

    The diffusion enhanced scheme by Bessemoulin-Chatard. Currently, the Pietra-Jüngel scheme is used for regularization of removable singularity.

    """

#todo_da: implemented diffusion enhanced scheme.
    function diffusionenhanced!(f, u, edge, data)
        tolRegularisation = 1.0e-13;

        uk  = viewK(edge, u)
        ul  = viewL(edge, u)

        ipsi = data.numberOfSpecies
        ireg = edge.region

        UT   = (kB * data.temperature ) / q
        dpsi = ul[ipsi]- uk[ipsi]
        f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

        for icc = 1:data.numberOfSpecies-1

            j0   = - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * UT * data.densityOfStates[ireg,icc]

            # etak = data.chargeNumbers[icc] / UT * ( uk[icc]-uk[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
            # etal = data.chargeNumbers[icc] / UT * ( ul[icc]-ul[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
            etak = etaFunction(uk,edge,data,icc,ipsi)
            etal = etaFunction(ul,edge,data,icc,ipsi)

            if abs( (etal-etak)/ (etak+etal)) > tolRegularisation
                g = (etal - etak ) / ( log(data.F(etal)) - log(data.F(etak)) )
            else # regularization idea coming from Pietra-Jüngel scheme
                gk = exp(etak)/data.F(etak)
                gl = exp(etal)/data.F(etal)
                g  = 0.5 * ( gk + gl )
            end

        nodel = edge.node[2]; nodek = edge.node[1]
        bandEdgeDifference = data.bandEdgeEnergyNode[nodel, icc] - data.bandEdgeEnergyNode[nodek, icc]

        bp, bm = fbernoulli_pm( data.chargeNumbers[icc] * (dpsi - bandEdgeDifference/q)/ (UT * g) )
        f[icc] = data.chargeNumbers[icc] * j0 * g * (  bm * data.F(etal) - bp * data.F(etak))
        end

    end

    """
    $(SIGNATURES)

    The Koprucki-Gärtner scheme. This scheme is calculated by solving a fixed point equation which arises when considering the generalized Scharfetter-Gummel scheme in case of Blakemore statistics.

    """
    # todo_da: is not working yet!! (solved problem!)
    function kopruckigaertner!(f, u, edge, data)
        gamma = 0.27        # from Blakemore distribution
        max_iteration = 200 # for Newton solver
        it = 0              # number of iterations (newton)
        damp = 0.1          # damping factor

        uk  = viewK(edge, u)
        ul  = viewL(edge, u)

        ipsi = data.numberOfSpecies
        ireg = edge.region

        UT   = (kB * data.temperature ) / q
        dpsi = ul[ipsi]- uk[ipsi]
        f[ipsi]  =  - data.dielectricConstant[ireg] * ε0 * dpsi

        for icc = 1:data.numberOfSpecies-1

            j0   = - data.chargeNumbers[icc] * q * data.mobility[ireg,icc] * UT * data.densityOfStates[ireg,icc]

            # etak = data.chargeNumbers[icc] / UT * ( uk[icc]-uk[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
            # etal = data.chargeNumbers[icc] / UT * ( ul[icc]-ul[ipsi] + data.bandEdgeEnergy[ireg,icc] / q)
            etak = etaFunction(uk,edge,data,icc,ipsi)
            etal = etaFunction(ul,edge,data,icc,ipsi)

            # use Sedan flux as starting guess
            Q = data.chargeNumbers[icc] * dpsi/ UT + (etal - etak) - log(data.F(etal)) + log(data.F(etak))
            bp, bm = fbernoulli_pm(Q)
            jInitial =  ( bm * data.F(etal)  - bp * data.F(etak))

#todo_da: nachvollziehen warum - gamma*j und nicht +
            implicitEquation(j::Real) =  (fbernoulli_pm(data.chargeNumbers[icc] * (dpsi / UT) - gamma*j )[2] * exp(etal) - fbernoulli_pm(data.chargeNumbers[icc] * (dpsi/ UT) - gamma*j )[1] * exp(etak)) - j

                delta = 1.0e-18 + 1.0e-14 * abs(value(jInitial))
                while (it < max_iteration)
                    Fval  = implicitEquation(jInitial)
                    dFval = ForwardDiff.derivative(implicitEquation,jInitial)
                    if isnan( value(dFval) ) || value( abs(dFval) ) < delta
                        @show value(j), value(Fval), value(dFval)
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

    Compute the electro-neutral solution, assuming the Boltzmann approximation. It is obtained by setting the left-hand side in
    the Poisson equation equal to zero and solving for \\psi.

    """
    #todo_da: Here: Need a grid from ExtendableGrids.jl
    #function electroNeutralSolutionBoltzmann(grid::VoronoiFVM.AbstractGrid,data::DDFermiData)
    # need to specify grid
    function electroNeutralSolutionBoltzmann(grid,data::DDFermiData)

        if data.numberOfSpecies-1 != 2
            error("The electroneutral solution is only implemented for two species!")
        end

        # region independent parameters
        iphin = 1;
        iphip = 2;
        UT    = (kB * data.temperature ) / q

        # initialize zero vector
    #todo_da:
    #    psi0 = zeros(length(grid.coord))
    coord        = grid[Coordinates]
    bfacenodes   = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    psi0         = zeros(length(coord))

        # boundary values
        #todo_da:
    #    for i=1:length(grid.bfacenodes)
        for i=1:length(bfacenodes)
            # boundary index
            # from old code: ibreg = grid.bfaceregions[i]
            ibreg = bfaceregions[i]
            # boundary region specific data
            Ec = data.bBandEdgeEnergy[ibreg,iphin]
            Ev = data.bBandEdgeEnergy[ibreg,iphip]
            Nc = data.bDensityOfStates[ibreg,iphin]
            Nv = data.bDensityOfStates[ibreg,iphip]
            Ni = sqrt( Nc*Nv*exp(-(Ec-Ev)/(kB*data.temperature)) )
            C  = -data.chargeNumbers[iphin] * data.bDoping[ibreg,iphin] -data.chargeNumbers[iphip] * data.bDoping[ibreg,iphip]
    #todo_da: Vorschlag: das minus aus doping nach vorne ziehen, damit klarer wird, dass es insgesamt subtrahiert wird (?)

            # set boundary values for electroneutral potential
psi0[bfacenodes[i]] = (Ec+Ev)/(2q) - 0.5*UT*log(Nc/Nv) + UT*asinh(C/(2*Ni))
            #todo_da: old code: psi0[grid.bfacenodes[i]] = (Ec+Ev)/(2q) - 0.5*UT*log(Nc/Nv) + UT*asinh(C/(2*Ni))
        end

        # interior values
        #todo_da
        cellregions = grid[CellRegions]
        #for i=1:length(grid.cellregions)-1
        for i=1:length(cellregions)-1

            # interior index
            ireg      = cellregions[i]
            ireg_next = cellregions[i+1]
    #todo_da: Das ist mir nicht ganz klar. Vermutung: Muss aufgteilt werden wegen Grenzschichten wahrscheinlich, weil da möglicherweise unterschiedliche Werte für physikalische Parameter. Aber warum das arithmetische Mittel dann von beiden? (Im Handbook steht nichts dazu)

            # interior region specific data
            Ec = (data.bandEdgeEnergy[ireg,iphin]  + data.bandEdgeEnergy[ireg_next,iphin] ) / 2 + data.bandEdgeEnergyNode[i,iphin]
            Ev = (data.bandEdgeEnergy[ireg,iphip]  + data.bandEdgeEnergy[ireg_next,iphip] ) / 2 + data.bandEdgeEnergyNode[i,iphip]
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


    """

    $(SIGNATURES)

    Find the equilibrium solution for the electrostatic potential

    """
    function solveEquilibriumBoltzmann!(solution, initialGuess, data, grid, control, dense)

        if !(data.F == DDFermi.Boltzmann) # if F != Boltzmann, find equilibrium solution for Boltzmann
#todo_da: add
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
        # todo_da: change in println and not prinln
            println("*** We compute with Boltzmann statistics anyway. ")

        end

    end
