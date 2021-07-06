##########################################################
##########################################################

"""
$(TYPEDEF)

A struct holding the physical region dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct ChargeTransportParams

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    """
    number of nodes used for the disretization of the domain ``\\mathbf{\\Omega}``.
    """
    numberOfNodes                ::  Int64

    """
    number of regions ``\\mathbf{\\Omega}_k`` within the domain ``\\mathbf{\\Omega}``.
    """
    numberOfRegions              ::  Int64

    """
    number of boundary regions ``(\\partial \\mathbf{\\Omega})_k`` such that 
    `` \\partial \\mathbf{\\Omega} = \\cup_k (\\partial \\mathbf{\\Omega})_k``.
    """
    numberOfBoundaryRegions      ::  Int64

    """
    number of moving charge carriers.
    """ 
    numberOfCarriers             ::  Int64


    """
    number of present inner interface carriers.
    """ 
    numberOfInterfaceCarriers    :: Int64


    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    """
    A given constant temperature.
    """
    temperature                  ::  Float64

    """
    The thermal voltage, which reads  ``U_T = k_B T / q``.
    """
    UT                           ::  Float64

    """
    A reference energy, which is only used for numerical computations.
    """
    # DA: I think that, when using Schottky contacts we need to set the reference
    # energy to the Fermi Level to get a physical reasonable electrostatic potential,
    # but I am not sure yet and did not test it well ...
    Eref                         ::  Float64

    """
    The parameter of the Blakemore statistics.
    """
    γ                            ::  Float64

    """
    Prefactor of electro-chemical reaction of internal boundary conditions.
    """
    r0                           ::  Float64


    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################
    """
    An array for the given applied voltages at the contacts.
    """
    contactVoltage               ::  Array{Float64,1}

    """
    An array for the given Fermi level at the contacts.
    """
    bFermiLevel                  ::  Array{Float64,1}


    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    """
    An array with the corresponding charge numbers
    ``z_\\alpha`` for all carriers ``\\alpha``.
    """
    chargeNumbers                ::  Array{Float64,1}
    

    ###############################################################
    ####    number of boundary regions x number of carriers    ####
    ###############################################################
    """
    A 2D array with the corresponding boundary band-edge energy values
    ``E_\\alpha`` for each carrier ``\\alpha``.
    """
    bBandEdgeEnergy              ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding boundary effective density of states values
    ``N_\\alpha`` for each carrier ``\\alpha``.
    """
    bDensityOfStates             ::  Array{Float64,2}

    """
    A 2D array with the corresponding boundary doping values for each carrier ``\\alpha``.
    """
    bDoping                      ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding boundary velocity values for each carrier ``\\alpha``,
    when assuming Schottky contacts.
    """
    bVelocity                    ::  Array{Float64,2}


    ###############################################################
    ####   number of bregions x 2 (for electrons and holes!)   ####
    ############################################################### 
    """
    A 2D array with the corresponding recombination surface boundary velocity values
    for electrons and holes.
    """
    recombinationSRHvelocity     ::  Array{Float64,2}
    """
    A 2D array with the corresponding recombination surface boundary density values
    for electrons and holes.
    """
    brecombinationSRHTrapDensity ::  Array{Float64,2}
 

    ###############################################################
    ####        number of regions x number of carriers         ####
    ###############################################################
    """
    A 2D array with the corresponding doping values for each carrier ``\\alpha`` on each region.
    """
    doping                       ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    densityOfStates              ::  Array{Float64,2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    bandEdgeEnergy               ::  Array{Float64,2}

    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha``
    for each carrier ``\\alpha`` on each region.
    """
    mobility                     ::  Array{Float64,2}


    ###############################################################
    #### number of regions x 2 (for electrons and holes only!) ####
    ############################################################### 
    """
    A 2D array with the corresponding SRH lifetimes ``\\tau_n, \\tau_p`` for electrons and holes.
    """
    recombinationSRHLifetime     ::  Array{Float64,2}

    """
    A 2D array with the corresponding SRH trap densities ``n_{\\tau}, p_{\\tau}`` for electrons and holes.
    """
    recombinationSRHTrapDensity  ::  Array{Float64,2}

    """
    A 2D array with the corresponding Auger coefficients for electrons and holes.
    """
    recombinationAuger           ::  Array{Float64,2}


    ###############################################################
    ####                   number of regions                   ####
    ############################################################### 
    """
    A region dependent dielectric constant.
    """
    dielectricConstant           ::  Array{Float64,1}

    """
    A region dependent array for the prefactor in the generation process which is the
    incident photon flux.
    """
    generationIncidentPhotonFlux ::  Array{Float64,1}
    """
    A region dependent array for an uniform generation rate.
    """
    generationUniform            ::  Array{Float64,1}

    """
    A region dependent array for the absorption coefficient in the generation process.
    """
    generationAbsorption         ::  Array{Float64,1}

    """
    A region dependent array for the radiative recombination rate.
    """
    recombinationRadiative       ::  Array{Float64,1}
    
    ###############################################################
    ChargeTransportParams() = new() # standard constructor

end


"""
$(TYPEDEF)

A struct holding the physical nodal, i.e. space, dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct ChargeTransportParamsNodal
    
    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    """
    A node dependent dielectric constant.
    """
    dielectricConstant           ::  Array{Float64,1}
 
 
    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    mobility                     ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding doping values for each carrier ``\\alpha`` on each node.
    """
    doping                       ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    densityOfStates              ::  Array{Float64,2}
 
    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha`` for each carrier ``\\alpha`` on each node.
    """
    bandEdgeEnergy               ::  Array{Float64,2}

    ###############################################################
    ChargeTransportParamsNodal() = new()

end

"""
$(TYPEDEF)

A struct holding all data information including model and numerics information,
but also all physical parameters for a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct ChargeTransportData

    ###############################################################
    ####                   model information                   ####
    ###############################################################
    """
    An array with the corresponding distribution function ``\\mathcal{F}_\\alpha`` for all carriers ``\\alpha``.
    """
    F                            ::  Array{Function,1}

    """
    An array of DataTypes with the type of boundary model for each boundary (interior and exterior).
    """
    boundary_type                ::  Array{DataType, 1}  

    """
    A DataType for the bulk recombination model.
    """
    bulk_recombination_model     ::  DataType   


    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    """
    A DataType for the flux discretization method.
    """
    flux_approximation           ::  DataType

    """
    A DataType for equilibrium or out of equilibrium calculations.
    """
    calculation_type             :: DataType

    """
    A DataType for transient or stationary calculations.
    """
    model_type                   :: DataType

    """
    A DataType for for generation model.
    """
    generation_model             :: DataType

    """
    An embedding parameter used to solve the nonlinear Poisson problem, which results
    in the case of thermodynamic equilibrium and electrocharge neutrality.
    """
    λ1                           ::  Float64

    """
    An embedding parameter for turning the generation rate ``G`` on.
    """
    λ2                           ::  Float64

    """
    An embedding parameter for electrochemical reaction.
    """
    λ3                           ::  Float64


    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    """
    A struct holding all region dependent parameter information. For more information see
    struct ChargeTransportParams.
    """
    params                       :: ChargeTransportParams

    """
    A struct holding all space dependent parameter information. For more information see
    struct ChargeTransportParamsNodal.
    """
    paramsnodal                  :: ChargeTransportParamsNodal

    ###############################################################
    ChargeTransportData() = new()

end


"""
$(TYPEDEF)

A struct holding all information necessary for a drift-diffusion type system.

$(TYPEDFIELDS)

"""
mutable struct ChargeTransportSystem

    """
    A struct holding all data information, see ChargeTransportData
    """
    data                         :: ChargeTransportData

    """
    A struct holding system information for the finite volume system.
    """
    fvmsys                       :: VoronoiFVM.AbstractSystem

    ###############################################################
    ChargeTransportSystem() = new()

end


##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Simplified constructor for ChargeTransportParams which only takes the grid
and the numberOfCarriers as argument.
    
"""
function ChargeTransportParams(grid, numberOfCarriers)

    numberOfNodes           = length(grid[Coordinates])
    numberOfRegions         = grid[NumCellRegions]
    numberOfBoundaryRegions = grid[NumBFaceRegions]
    ###############################################################

    params = ChargeTransportParams()

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    params.numberOfNodes                = numberOfNodes
    params.numberOfRegions              = numberOfRegions
    params.numberOfBoundaryRegions      = numberOfBoundaryRegions
    params.numberOfCarriers             = numberOfCarriers
    params.numberOfInterfaceCarriers    = 0

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    params.temperature                  = 300 * K                 
    params.UT                           = (kB * 300 * K ) / q     # thermal voltage
    params.Eref                         = 0.0                     # reference energy
    params.γ                            = 0.27                    # parameter for Blakemore statistics
    params.r0                           = 0.0                     # r0 prefactor electro-chemical reaction

    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################
    params.contactVoltage               = spzeros(Float64, numberOfBoundaryRegions)
    params.bFermiLevel                  = spzeros(Float64, numberOfBoundaryRegions)    # Fermi level at boundary
    
    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    params.chargeNumbers                = spzeros(Float64, numberOfCarriers) 

    ###############################################################
    ####    number of boundary regions x number of carriers    ####
    ###############################################################
    params.bBandEdgeEnergy              = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDensityOfStates             = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDoping                      = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bVelocity                    = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)  # velocity at boundary; SchottkyContacts

    ###############################################################
    ####   number of bregions x 2 (for electrons and holes!)   ####
    ############################################################### 
    params.brecombinationSRHTrapDensity = spzeros(Float64, 2, numberOfBoundaryRegions)                 # for surface reco
    params.recombinationSRHvelocity     = spzeros(Float64, 2, numberOfBoundaryRegions)                 # for surface reco

    ###############################################################
    ####        number of regions x number of carriers         ####
    ###############################################################
    params.doping                       = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.densityOfStates              = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.bandEdgeEnergy               = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.mobility                     = spzeros(Float64, numberOfCarriers, numberOfRegions)

    ###############################################################
    #### number of regions x 2 (for electrons and holes only!) ####
    ###############################################################
    params.recombinationSRHLifetime     = Array{Float64,2}(undef, 2, numberOfRegions)
    params.recombinationSRHTrapDensity  = Array{Float64,2}(undef, 2, numberOfRegions)
    params.recombinationAuger           = Array{Float64,2}(undef, 2, numberOfRegions)

    ###############################################################
    ####                   number of regions                   ####
    ############################################################### 
    params.dielectricConstant           = spzeros(Float64, numberOfRegions)
    params.generationUniform            = spzeros(Float64, numberOfRegions)
    params.generationIncidentPhotonFlux = spzeros(Float64, numberOfRegions)
    params.generationAbsorption         = spzeros(Float64, numberOfRegions)
    params.recombinationRadiative       = spzeros(Float64, numberOfRegions)

    ###############################################################
    return params

end

"""
$(TYPEDSIGNATURES)

Simplified constructor for ChargeTransportParamsNodal which only takes the grid
and the numberOfCarriers as argument.
    
"""
function ChargeTransportParamsNodal(grid, numberOfCarriers)

    numberOfNodes                       = length(grid[Coordinates])

    ###############################################################

    paramsnodal                         = ChargeTransportParamsNodal()

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    paramsnodal.dielectricConstant      = spzeros(Float64, numberOfNodes) 
 
 
    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    paramsnodal.mobility                = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.doping                  = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.densityOfStates         = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.bandEdgeEnergy          = spzeros(Float64, numberOfCarriers, numberOfNodes)

    ###############################################################
    return paramsnodal

end


"""
$(TYPEDSIGNATURES)

Simplified constructor for ChargeTransportData which only takes the grid
and the numberOfCarriers as argument. Here all necessary information
including the physical parameters are located.
    
"""
function ChargeTransportData(grid, numberOfCarriers)

    numberOfBoundaryRegions = grid[NumBFaceRegions]

    ###############################################################
    data = ChargeTransportData()

    ###############################################################
    ####                   model information                   ####
    ###############################################################
    #data.indexSet  = set_indices(grid, numberOfCarriers, interface_model_none) # without interface model 

    data.F                        = fill!(similar(Array{Function,1}(undef, numberOfCarriers),Function), Boltzmann)
    data.boundary_type            = Array{DataType,1}(undef, numberOfBoundaryRegions)
    data.bulk_recombination_model = bulk_recombination_none

    for ii in 1:numberOfBoundaryRegions # as default all boundaries are set to an ohmic contact model.
        data.boundary_type[ii] = ohmic_contact
    end

    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    data.flux_approximation       = ScharfetterGummel
    data.calculation_type         = inEquilibrium           # do performances inEquilibrium or outOfEquilibrium
    data.model_type               = model_stationary        # indicates if we need additional time dependent part
    data.generation_model         = generation_none         # generation model
    data.λ1                       = 0.0                     # λ1: embedding parameter for NLP
    data.λ2                       = 0.0                     # λ2: embedding parameter for G
    data.λ3                       = 0.0                     # λ3: embedding parameter for electro chemical reaction

    
    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    data.params                   = ChargeTransportParams(grid, numberOfCarriers)
    data.paramsnodal              = ChargeTransportParamsNodal(grid, numberOfCarriers)
 
    ###############################################################

    return data

end

"""
$(TYPEDSIGNATURES)

Simplified constructor for ChargeTransportSystem. This is the main
struct in which all information is stored and with which the calculations
are performed.
    
"""
function ChargeTransportSystem(grid, data ;unknown_storage)

    num_species  = data.params.numberOfCarriers + data.params.numberOfInterfaceCarriers + 1

    ctsys        = ChargeTransportSystem()

    ctsys.data   = data
    
    physics      = VoronoiFVM.Physics(data        = data,
                                      num_species = num_species,
                                      flux        = flux!,
                                      reaction    = reaction!,
                                      breaction   = breaction!,
                                      storage     = storage!,
                                      bstorage    = bstorage!
                                      )

    ctsys.fvmsys = VoronoiFVM.System(grid, physics, unknown_storage = unknown_storage)

    # for detection of number of species. In following versions we can simply delete num_species from physics initialization. 
    VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species) 

    return ctsys

end

function show_params(ctsys::ChargeTransportSystem)

    params = ctsys.data.params
    for name in fieldnames(typeof(params))[1:end] 
        @printf("%30s = ",name)
        println(getfield(params,name))
    end

end


function Base.show(io::IO, this::ChargeTransportParamsNodal)
    for name in fieldnames(typeof(this))[1:end] 
        @printf("%30s = ",name)
        println(io,getfield(this,name))
    end
end
###########################################################
###########################################################
"""
$(TYPEDSIGNATURES)

Functions which sets for given charge carrier at a given boundary
a given value.
    
"""
function set_ohmic_contact!(ctsys, icc, ibreg, contact_val)
 
    ctsys.fvmsys.boundary_factors[icc, ibreg] = VoronoiFVM.Dirichlet
    ctsys.fvmsys.boundary_values[icc, ibreg]  = contact_val
 
end

###########################################################
###########################################################
# Wrappers for methods of VoronoiFVM

ct_enable_species!(ctsys::ChargeTransportSystem, ispecies, regions)   = VoronoiFVM.enable_species!(ctsys.fvmsys, ispecies, regions)

ct_enable_boundary_species!(ctsys::ChargeTransportSystem, ispecies, regions)   = VoronoiFVM.enable_boundary_species!(ctsys.fvmsys, ispecies, regions)

ct_unknowns(ctsys::ChargeTransportSystem)                             = VoronoiFVM.unknowns(ctsys.fvmsys)

ct_solve!(solution, initialGuess, ctsys, ;control=control, tstep=Inf) = VoronoiFVM.solve!(solution, initialGuess, ctsys.fvmsys, control=control, tstep=tstep)

###########################################################
###########################################################
"""
$(TYPEDSIGNATURES)

Functions which sets for given charge carrier at a given boundary
a given value.
    
"""

function ct_equilibrium_solve!(ctsys::ChargeTransportSystem; control = VoronoiFVM.NewtonControl(), nonlinear_steps = 20.0)

    ctsys.data.calculation_type    = inEquilibrium

    # initialize solution and starting vectors
    initialGuess                   = ct_unknowns(ctsys)
    solution                       = ct_unknowns(ctsys)
    @views initialGuess           .= 0.0 

    # we slightly turn a linear Poisson problem to a nonlinear one with these variables.
    I      = collect(nonlinear_steps:-1:0.0)
    LAMBDA = 10 .^ (-I) 
    prepend!(LAMBDA, 0.0)

    for i in 1:length(LAMBDA)

        if control.verbose
            println("λ1 = $(LAMBDA[i])")
        end
        ctsys.fvmsys.physics.data.λ1 = LAMBDA[i] 
        try

            ct_solve!(solution, initialGuess, ctsys, control = control, tstep=Inf)

        catch
            if (control.handle_exceptions)
                error("try to adjust nonlinear_steps, currently set to $(nonlinear_steps) or adjust Newton control parameters.")
            end
        end
    
        initialGuess = solution
    end

    return solution

end
###########################################################
###########################################################
function set_indices!(grid, numberOfCarriers, ::Type{interface_model_none}) # we are in most classical setting

    indexSet = Dict()

    indexSet["iphin"] = 1
    indexSet["iphip"] = 2

    if numberOfCarriers == 3

        indexSet["iphia"] = 3

    end
    
    indexSet["ipsi"]  = numberOfCarriers + 1

    return indexSet
end

function set_indices!(grid, numberOfCarriers, ::Type{interface_model_ion_charge})
# DA: generalizing to arbitrary domains is a bit difficult since we do not know which ones 
# are the active perovskite layers. Thus, we additionally need them here as input argument?
    bcellregions = grid[NumBFaceRegions]
    indexSet = Dict()

    indexSet["iphin"] = 1
    indexSet["iphip"] = 2
    indexSet["iphia"] = 3


    indexSet["iphiaJunction"] = 4:5
    indexSet["ipsi"] = 3 + length( indexSet["iphiaJunction"] ) + 1
        
    return indexSet

end


function set_indices!(grid, numberOfCarriers, ::Type{interface_model_surface_recombination})

    cellregions = grid[NumCellRegions]
    indexSet = Dict()

    indexSet["iphin"] = collect(1:cellregions)
    indexSet["iphip"] = collect(cellregions+1:2*cellregions)

    if numberOfCarriers == 3
        indexSet["iphia"] = 2*cellregions+1
    elseif numberOfCarriers < 2 || numberOfCarriers > 3
        println("Case of more than three carriers not tested yet. Hence, only electrons and holes are assumed.")
    end

    indexSet["ipsi"] = ( numberOfCarriers+1 ) + 2 * (cellregions-1)
        
    return indexSet

end
###########################################################
###########################################################
"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities.

"""
function compute_densities!(u, data, node, region, icc::Int, ipsi::Int, in_region::Bool)

    params      = data.params
    paramsnodal = data.paramsnodal 

    if in_region == false
        (params.bDensityOfStates[icc, region] + paramsnodal.densityOfStates[icc, node] ) * data.F[icc](etaFunction(u, data, node, region, icc, ipsi, in_region::Bool))
    elseif in_region == true
        (params.densityOfStates[icc, region] + paramsnodal.densityOfStates[icc, node])* data.F[icc](etaFunction(u, data, node, region, icc, ipsi, in_region::Bool)) 
    end
        
end


"""

$(TYPEDSIGNATURES)

For given potentials in vector form, compute corresponding vectorized densities.
[Caution: this was not tested for multidimensions.]
"""
function compute_densities!(grid, data, sol)
    params       = data.params
    paramsnodal  = data.paramsnodal

    ipsi         = params.numberOfCarriers + 1 
    densities    = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
        
    bfacenodes   = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    cellRegions  = copy(grid[CellRegions])
    cellRegions  = push!(cellRegions, grid[CellRegions][end]) #  enlarge region by final cell
    
    if dim_space(grid) > 1
        println("compute_densities! is so far only tested in 1D")
    end
    
    for icc in 1:params.numberOfCarriers

        for node in 1:params.numberOfNodes
            in_region = true
            u         = sol[:, node]
            region    = cellRegions[node]

            if node in bfacenodes
                in_region = false
                indexNode = findall(x -> x == node, vec(bfacenodes))[1]  # we need to know which index the node has in bfacenodes
                region    = bfaceregions[indexNode]                      # since the corresponding region number is at the same index
            end

            densities[icc, node] = compute_densities!(u, data, node, region, icc, ipsi, in_region)
        end
    
    end
    
    return densities

end


"""

$(SIGNATURES)

For given solution in vector form, compute corresponding vectorized band-edge energies and Fermi level.
[Caution: this was not tested for multidimensions.]
"""
function compute_energies!(grid, data, sol)

    params       = data.params
    paramsnodal  = data.paramsnodal
    
    ipsi         = params.numberOfCarriers + 1
    energies     = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
    fermiLevel   = Array{Real,2}(undef, data.numberOfCarriers, size(sol, 2))
    
    cellregions  = grid[CellRegions]
    cellregions  = push!(cellregions, cellregions[end])
    
    for icc in 1:params.numberOfCarriers

        for inode in 1:params.numberOfNodes
             E                      = params.bandEdgeEnergy[icc, cellregions[inode]] + paramsnodal.bandEdgeEnergy[icc, inode]
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
function electroNeutralSolution!(grid, data; Newton=false)

    params          = data.params
    paramsnodal     = data.paramsnodal

    if params.numberOfCarriers > 2
        error("this method is currently only working for electrons and holes")
    end

    solution        = zeros(length(grid[Coordinates]))
    iccVector       = collect(1:params.numberOfCarriers)
    zVector         = params.chargeNumbers[iccVector]
    FVector         = data.F[iccVector]
    regionsAllCells = copy(grid[CellRegions])
    regionsAllCells = push!(regionsAllCells, grid[CellRegions][end]) #  enlarge region by final cell
    phi             = 0.0                                            # in equilibrium set to 0
    psi0_initial    = 0.5

    for index = 1:length(regionsAllCells) - 1
        
        ireg          = regionsAllCells[index]
        zVector       = params.chargeNumbers[iccVector]
        FVector       = data.F[iccVector]
        regionsOfCell = regionsAllCells[grid[CellNodes][:,index]]   # all regions of nodes belonging to cell for given index

        # average following quantities if needed among all regions
        EVector = Float64[]; CVector = Float64[]; NVector = Float64[]

        for icc = 1:params.numberOfCarriers
            push!(EVector, sum(params.bandEdgeEnergy[icc, regionsOfCell])  / length(regionsOfCell) + paramsnodal.bandEdgeEnergy[icc,index])
            push!(CVector, sum(params.doping[icc, regionsOfCell])          / length(regionsOfCell) + paramsnodal.doping[icc, index])
            push!(NVector, sum(params.densityOfStates[icc, regionsOfCell]) / length(regionsOfCell) + paramsnodal.densityOfStates[icc, index])
        end
        # rhs of Poisson's equation as anonymous function depending on psi0
        f = psi0 -> chargeDensity(psi0, phi, params.UT, EVector, zVector, CVector, NVector, FVector)

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
