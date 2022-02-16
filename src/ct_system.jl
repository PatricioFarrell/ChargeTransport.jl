##########################################################
##########################################################

"""
$(TYPEDEF)

A struct holding all information necessary for building bulk recombination.
With help of this constructor we can read out the indices the user chooses for
electron and hole quasi Fermi potentials.

$(TYPEDFIELDS)

"""
mutable struct BulkRecombination

    """
    Index for data construction of quasi Fermi potential of electrons
    """
	iphin                 ::  Int64

    """
    Index for data construction of quasi Fermi potential of holes
    """
    iphip                 ::  Int64

    """
    Boolean for present Auger recombination in bulk
    """
    bulk_recomb_Auger     ::  Bool

    """
    Boolean for present radiative recombination in bulk
    """
    bulk_recomb_radiative ::  Bool

    """
    DataType for present SRH recombination in bulk
    """
    bulk_recomb_SRH       ::  Union{SRHWithoutTraps, SRHWithTraps}

    BulkRecombination() = new()

end

"""
$(TYPEDEF)

A struct holding all information necessary for Interface species.
With help of this constructor we can read out the indices and boundaries
the user chooses for interface species.

$(TYPEDFIELDS)

"""

"""
$(SIGNATURES)

Corresponding constructor for the bulk recombination model.
"""
function set_bulk_recombination(;iphin = 1, iphip = 2,
                                  bulk_recomb_Auger = true,
                                  bulk_recomb_radiative = true,
                                  bulk_recomb_SRH = true)

    bulk_recombination = BulkRecombination()

    bulk_recombination.iphin                 = iphin
    bulk_recombination.iphip                 = iphip

    bulk_recombination.bulk_recomb_Auger     = bulk_recomb_Auger
    bulk_recombination.bulk_recomb_radiative = bulk_recomb_radiative

    if bulk_recomb_SRH == true
        bulk_recombination.bulk_recomb_SRH   = SRHStationary
    else
        bulk_recombination.bulk_recomb_SRH   = SRHOff
    end

    return bulk_recombination

end

"""
$(TYPEDEF)

A struct holding all information necessary for enabling traps
in the SRH recombination.
With help of this constructor we can read out the index the user chooses for
trap quasi Fermi potentials and the respective regions in which they are
defined.

$(TYPEDFIELDS)

"""
mutable struct Traps

    """
    Array with the index of traps.
    """
    traps       ::  Int64

    """
    Corresponding regions where traps are assumed to be present.
    """
    regions     ::  Array{Int64, 1}

    Traps() = new()

end

"""
$(TYPEDEF)

Corresponding constructor for the present trap density and the respective regions.
"""
function enable_traps!(;data = data, traps = 3, regions = [1, 2, 3])

    enable_traps                            = Traps()

    enable_traps.traps                      = traps
    enable_traps.regions                    = regions

    if data.model_type == Transient
        data.bulk_recombination.bulk_recomb_SRH = SRHTrapsTransient
    else
        data.bulk_recombination.bulk_recomb_SRH = SRHStationary
    end

    data.enable_traps = enable_traps

end

###########################################################
###########################################################

"""
$(TYPEDEF)

A struct holding all information necessary on the ionic charge carriers.
With help of this constructor we can read out the indices the user chooses for
ionic charge carrier quasi Fermi potentials and the respective regions in which they are
defined. Note that it is possible to use ions as well as ion vacancies.

$(TYPEDFIELDS)

"""
mutable struct IonicChargeCarriers

    """
    Array with the indices of ionic charge carriers.
    """
    ionic_carriers       ::  Array{Int64, 1}

    """
    Corresponding regions where ionic charge carriers are assumed to be present.
    """
    regions              ::  Array{Int64, 1}

    IonicChargeCarriers() = new()

end


"""
Corresponding constructor for the present ionic charge carriers and the respective regions.
"""
function enable_ionic_carriers(;ionic_carriers = [3], regions = [2])

    enable_ions = IonicChargeCarriers()

    enable_ions.ionic_carriers   = ionic_carriers
    enable_ions.regions          = regions

    return enable_ions

end


###########################################################
###########################################################

mutable struct InterfaceCarriers

    """
    index for data construction of Interface species
    """
	index                 ::  Array{Int64,1} # here the id's of the AbstractQuantities or the integer indices are parsed.

    """
    boundary region number
    """
    boundary_region       ::  Int64

    InterfaceCarriers() = new()

end

"""
$(SIGNATURES)

Corresponding constructor for the interface species.
"""
function enable_interface_carrier!(data,; species::Array{Int64, 1}, boundary_region::Int64)

    interfaceCarriers = InterfaceCarriers()

    interfaceCarriers.index           = species
    interfaceCarriers.boundary_region = boundary_region

    data.interfaceCarriers            = interfaceCarriers

end

###########################################################
###########################################################

"""
$(TYPEDEF)

A struct holding the physical region dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct Params

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    """
    Number of nodes used for the disretization of the domain ``\\mathbf{\\Omega}``.
    """
    numberOfNodes                ::  Int64

    """
    Number of subregions ``\\mathbf{\\Omega}_k`` within the domain ``\\mathbf{\\Omega}``.
    """
    numberOfRegions              ::  Int64

    """
    Number of boundary regions ``(\\partial \\mathbf{\\Omega})_k`` such that
    `` \\partial \\mathbf{\\Omega} = \\cup_k (\\partial \\mathbf{\\Omega})_k``.
    Note that here are interior and outerior boundaries calculated.
    """
    numberOfBoundaryRegions      ::  Int64

    """
    Number of moving charge carriers.
    """
    numberOfCarriers             ::  Int64

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
    An array for the given Schottky barriers at present Schotkky contacts.
    """
    SchottkyBarrier              ::  Array{Float64,1}

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
    An array with the corresponding boundary band-edge energy values
    ``E_\\alpha`` in each region for each carrier ``\\alpha``.
    """
    bBandEdgeEnergy              ::  Array{Float64,2}

    """
    An array with the corresponding boundary effective density of states values
    ``N_\\alpha`` for each carrier ``\\alpha``.
    """
    bDensityOfStates             ::  Array{Float64,2}


    """
    A 2D array with the corresponding boundary mobility values `` \\mu_\\alpha``
    in each boundary region for each carrier ``\\alpha``.
    """
    bMobility                    ::  Array{Float64,2}

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
    bRecombinationSRHTrapDensity ::  Array{Float64,2}


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
    A 2D array with the corresponding time-independent SRH trap densities
    ``n_{\\tau}, p_{\\tau}`` for electrons and holes.
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
    Params() = new() # standard constructor

end


"""
$(TYPEDEF)

A struct holding the physical nodal, i.e. space-dependent parameters for
a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct ParamsNodal

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    """
    A node dependent dielectric constant.
    """
    dielectricConstant           ::  Array{Float64,1}

    """
    A 1D array with the corresponding doping values on each node.
    """
    doping                       ::  Array{Float64,1}
    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    """
    A 2D array with the corresponding mobility values ``\\mu_\\alpha`` for each carrier
    ``\\alpha`` on each node.
    """
    mobility                     ::  Array{Float64,2}

    """
    A 2D array with the corresponding effective density of states values ``N_\\alpha`` for
    each carrier ``\\alpha`` on each node.
    """
    densityOfStates              ::  Array{Float64,2}

    """
    A 2D array with the corresponding band-edge energy values ``E_\\alpha`` for each carrier
    ``\\alpha`` on each node.
    """
    bandEdgeEnergy               ::  Array{Float64,2}

    ###############################################################
    ParamsNodal() = new()

end

"""
$(TYPEDEF)

A struct holding all data information including model and numerics information,
but also all physical parameters for a drift-diffusion simulation of a semiconductor device.

$(TYPEDFIELDS)

"""
mutable struct Data{TFuncs<:Function}

    # DA: this one stands for interface thickness and is introduced for test set-up
    d                            :: Float64
    ###############################################################
    ####                   model information                   ####
    ###############################################################
    """
    An array with the corresponding distribution function ``\\mathcal{F}_\\alpha`` for all
    carriers ``\\alpha``.
    """
    F                            ::  Array{TFuncs,1}

    """
    An array of DataTypes with the type of boundary model for each boundary (interior and exterior).
    """
    boundary_type                ::  Array{BoundaryModel, 1}

    """
    A struct containing information concerning the bulk recombination model.
    """
    bulk_recombination           ::  BulkRecombination

    """
    A struct which contains information on the regions, where ionic charge carriers
    (ions and/or ion vacancies) are present.
    """
    enable_ionic_carriers        ::  IonicChargeCarriers

    """
    An AbstractVector which contains information on present SRH traps.
    """
    enable_traps                 ::  Traps

    """
    A struct which contains information on present interface species.
    """
    interfaceCarriers            :: InterfaceCarriers

    """
    DataType which stores information about which inner interface model is chosen by user.
    This quantity cannot be seen by the user and is needed for the core of package.
    """
    inner_interface_model        ::  InterfaceModel

    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    """
    A DataType for the flux discretization method.
    """
    flux_approximation           ::  FluxApproximations

    """
    A DataType for equilibrium or out of equilibrium calculations.
    """
    calculation_type             ::  CalculationType

    """
    A DataType for transient or stationary calculations.
    """
    model_type                   ::  ModelType

    """
    A DataType for for generation model.
    """
    generation_model             ::  GenerationModel

    """
    An embedding parameter used to solve the nonlinear Poisson problem, where for
    λ1 = 0 the right hand-side is set to zero whereas for
    for λ1 = 1 we have a full space charge density.
    """
    λ1                           ::  Float64

    """
    An embedding parameter for the generation rate ``G``.
    """
    λ2                           ::  Float64

    """
    An embedding parameter for the electrochemical reaction.
    """
    λ3                           ::  Float64

    ###############################################################
    ####             Templates for DOS and BEE                 ####
    ###############################################################

    """
    Within this template informations concerning the band-edge energy
    of each carrier is stored locally which saves allocations.
    We have two of such templates due to the two point flux approximation schemes.
    """
    tempBEE1                     ::  Array{Float64, 1}

    """
    See the description of tempBEE1
    """
    tempBEE2                     ::  Array{Float64, 1}

    """
    Within this template informations concerning the effective DOS
    of each carrier is stored locally which saves allocations.
    We have two of such templates due to the two point flux approximation schemes.
    """
    tempDOS1                     ::  Array{Float64, 1}

    """
    See the desciption of tempDOS2
    """
    tempDOS2                     ::  Array{Float64, 1}

    ###############################################################
    ####        Quantities (for discontinuous solving)         ####
    ###############################################################

    """
    An array which contains information on whether charge carriers
    are continuous or discontinuous.
    This is needed for building the AbstractQuantities which handle the
    indices of charge carriers on different regions.
    """
    isContinuous                 :: Array{Bool, 1}


    """
    This list stores all charge carriers. Based on
    the user choice we have with this new type the opportunity to
    simulate discontinuous unknowns.
    """
    chargeCarrierList            :: Array{QType, 1}

    """
    This list stores all interface charge carriers.
    Here, we can have a vector holding all AbstractQuantities
    or a vector holding an integer array depending on the input.
    """
    interfaceCarrierList         :: Union{Array{VoronoiFVM.AbstractQuantity,1}, Array{Int64, 1}}

    """
    This variable stores the index of the electric potential. Based on
    the user choice we have with this new type the opportunity to
    simulate discontinuous unknowns.
    """
    index_psi                    :: QType

    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    """
    A struct holding all region dependent parameter information. For more information see
    struct Params.
    """
    params                       :: Params

    """
    A struct holding all space dependent parameter information. For more information see
    struct ParamsNodal.
    """
    paramsnodal                  :: ParamsNodal

    ###############################################################
    Data{TFuncs}() where {TFuncs} = new()

end



"""
$(TYPEDEF)

A struct holding all information necessary for a drift-diffusion type system.

$(TYPEDFIELDS)

"""
mutable struct System

    """
    A struct holding all data information, see Data
    """
    data                         :: Data

    """
    A struct holding system information for the finite volume system.
    """
    fvmsys                       :: VoronoiFVM.AbstractSystem

    ###############################################################
    System() = new()

end


##########################################################
##########################################################

"""
$(TYPEDSIGNATURES)

Simplified constructor for Params which only takes the grid
and the numberOfCarriers as argument.

"""
function Params(grid, numberOfCarriers)

    numberOfNodes           = size(grid[Coordinates], 2)
    numberOfRegions         = grid[NumCellRegions]
    numberOfBoundaryRegions = grid[NumBFaceRegions]
    ###############################################################

    params = Params()

    ###############################################################
    ####                   integer numbers                     ####
    ###############################################################
    params.numberOfNodes                = numberOfNodes
    params.numberOfRegions              = numberOfRegions
    params.numberOfBoundaryRegions      = numberOfBoundaryRegions
    params.numberOfCarriers             = numberOfCarriers

    ###############################################################
    ####                     real numbers                      ####
    ###############################################################
    params.temperature                  = 300 * K
    params.UT                           = (kB * 300 * K ) / q # thermal voltage
    params.γ                            = 0.27                # parameter for Blakemore statistics
    params.r0                           = 0.0                 # r0 prefactor electro-chemical reaction

    ###############################################################
    ####              number of boundary regions               ####
    ###############################################################
    params.SchottkyBarrier              = spzeros(Float64, numberOfBoundaryRegions)

    ###############################################################
    ####                  number of carriers                   ####
    ###############################################################
    params.chargeNumbers                = spzeros(Float64, numberOfCarriers)

    ###############################################################
    ####    number of boundary regions x number of carriers    ####
    ###############################################################
    params.bBandEdgeEnergy              = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDensityOfStates             = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bMobility                    = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bDoping                      = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.bVelocity                    = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)

    ###############################################################
    ####   number of bregions x 2 (for electrons and holes!)   ####
    ###############################################################
    params.bRecombinationSRHTrapDensity = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)
    params.recombinationSRHvelocity     = spzeros(Float64, numberOfCarriers, numberOfBoundaryRegions)

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
    params.recombinationSRHLifetime     = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationSRHTrapDensity  = spzeros(Float64, numberOfCarriers, numberOfRegions)
    params.recombinationAuger           = spzeros(Float64, numberOfCarriers, numberOfRegions)

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

Simplified constructor for ParamsNodal which only takes the grid
and the numberOfCarriers as argument.

"""
function ParamsNodal(grid, numberOfCarriers)

    numberOfNodes                       = length(grid[Coordinates])

    ###############################################################

    paramsnodal                         = ParamsNodal()

    ###############################################################
    ####                    number of nodes                    ####
    ###############################################################
    paramsnodal.dielectricConstant      = spzeros(Float64, numberOfNodes)
    paramsnodal.doping                  = spzeros(Float64, numberOfNodes)


    ###############################################################
    ####          number of nodes x number of carriers         ####
    ###############################################################
    paramsnodal.mobility                = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.densityOfStates         = spzeros(Float64, numberOfCarriers, numberOfNodes)
    paramsnodal.bandEdgeEnergy          = spzeros(Float64, numberOfCarriers, numberOfNodes)

    ###############################################################
    return paramsnodal

end


"""
$(TYPEDSIGNATURES)

Simplified constructor for Data which only takes the grid
and the numberOfCarriers as argument. Here, all necessary information
including the physical parameters, but also some numerical information
are located.

"""
function Data(grid, numberOfCarriers; statfunctions::Type{TFuncs}=StandardFuncSet) where TFuncs

    numberOfBoundaryRegions = grid[NumBFaceRegions]

    ###############################################################
    data = Data{TFuncs}()

    ###############################################################
    ####                   model information                   ####
    ###############################################################

    data.F                      = TFuncs[ Boltzmann for i=1:numberOfCarriers]
    data.boundary_type          = BoundaryModel[ InterfaceModelNone for i = 1:numberOfBoundaryRegions]

    # bulk_recombination is a struct holding the input information
    data.bulk_recombination     = set_bulk_recombination(iphin = 1, iphip = 2, bulk_recomb_Auger = false, bulk_recomb_radiative = false,
                                                         bulk_recomb_SRH = false)

    # enable_ionic_carriers is a struct holding the input information
    data.enable_ionic_carriers  = enable_ionic_carriers(ionic_carriers = [3], regions = [2])

    data.enable_traps           = Traps()

    data.inner_interface_model  = InterfaceModelNone

    ###############################################################
    ####                 Numerics information                  ####
    ###############################################################
    data.flux_approximation     = ScharfetterGummel
    data.calculation_type       = InEquilibrium       # do performances InEquilibrium or OutOfEquilibrium
    data.model_type             = Stationary          # indicates if we need additional time dependent part
    data.generation_model       = GenerationNone      # generation model
    data.λ1                     = 1.0                 # λ1: embedding parameter for NLP
    data.λ2                     = 1.0                 # λ2: embedding parameter for G
    data.λ3                     = 1.0                 # λ3: embedding parameter for electro chemical reaction

    ###############################################################
    ####             Templates for DOS and BEE                 ####
    ###############################################################

    data.tempBEE1               = spzeros(Float64, numberOfCarriers)
    data.tempBEE2               = spzeros(Float64, numberOfCarriers)
    data.tempDOS1               = spzeros(Float64, numberOfCarriers)
    data.tempDOS2               = spzeros(Float64, numberOfCarriers)

    ###############################################################
    ####        Quantities (for discontinuous solving)         ####
    ###############################################################
    data.isContinuous           = Bool[ true for i = 1:numberOfCarriers]
    # default values for most simple case
    data.chargeCarrierList      = collect(1:numberOfCarriers)

    data.index_psi              = numberOfCarriers + 1

    ###############################################################
    ####          Physical parameters as own structs           ####
    ###############################################################
    data.params                 = Params(grid, numberOfCarriers)
    data.paramsnodal            = ParamsNodal(grid, numberOfCarriers)

    ###############################################################

    return data

end

###########################################################
###########################################################

"""
$(SIGNATURES)

System constructor which builds all necessary information needed based
on the input parameters with special regard to additional interface models.
This is the main struct in which all information on the input data, but also
on the solving system, with which the calculations are performed,
are stored.

"""
function System(grid, data ;unknown_storage)

    ctsys                 = System()

    interface_model       = inner_interface_model(data)
    # here at this point, based on the interface model, we choose a system based on normal
    # integer indexing or quantity indexing.
    ctsys                 = build_system(grid, data, unknown_storage, interface_model)

    return ctsys

end


# The core of the system constructor. Here, the system for no additional interface model is build.
function build_system(grid, data, unknown_storage, ::Type{InterfaceModelNone})

    num_species_sys        = data.params.numberOfCarriers + 1

    ctsys                  = System()

    # save this information such that there is no need to calculate it again for boundary conditions
    data.inner_interface_model = InterfaceModelNone

    # Here, in this case for the loops within physics methods we set the chargeCarrierList to normal indexing.
    data.chargeCarrierList = collect(1:data.params.numberOfCarriers)
    # DA: caution with the interface_model with ionic interface charges (in future versions,
    # we will work with VoronoiFVM.InterfaceQuantites)

    # put Auger and radiative on or off
    if data.bulk_recombination.bulk_recomb_Auger == false
        data.params.recombinationAuger .= 0.0
    end

    if data.bulk_recombination.bulk_recomb_radiative == false
        data.params.recombinationRadiative .= 0.0
    end

    data.index_psi         = num_species_sys

    ctsys.data             = data

    physics      = VoronoiFVM.Physics(data        = data,
                                      flux        = flux!,
                                      reaction    = reaction!,
                                      breaction   = breaction!,
                                      storage     = storage!,
                                      bstorage    = bstorage!,
                                      bflux       = bflux!
                                      )

    ctsys.fvmsys = VoronoiFVM.System(grid, physics, unknown_storage = unknown_storage)

    for icc ∈ data.chargeCarrierList # first simply define all charge carriers on whole domain
        enable_species!(ctsys.fvmsys, icc, 1:data.params.numberOfRegions)
    end

    if data.params.numberOfCarriers > 2 # when ionic vacancies are present

        # get the ion vacancy indices by user input which where parsed into the struct
        # data.enable_ionic_carriers
        ionic_carriers = data.enable_ionic_carriers.ionic_carriers

        for icc ∈ ionic_carriers # Then, ion vacancies only present in user defined layers, i.e. adjust previous specification
            enable_species!(ctsys.fvmsys, icc, data.enable_ionic_carriers.regions)
        end
    end

    enable_species!(ctsys.fvmsys, data.index_psi, 1:data.params.numberOfRegions) # ipsi defined on whole domain

    # for detection of number of species
    VoronoiFVM.increase_num_species!(ctsys.fvmsys, num_species_sys)

    return ctsys

end


# The core of the new system constructor. Here, the system for discontinuous quantities
# is build.
function build_system(grid, data, unknown_storage, ::Type{InterfaceModelDiscontqF})

    ctsys        = System()
    fvmsys       = VoronoiFVM.System(grid, unknown_storage=unknown_storage)

    # save this information such that there is no need to calculate it again for boundary conditions
    data.inner_interface_model = InterfaceModelDiscontqF

    # DA: Not working like that anymore ......
    # if data.params.numberOfCarriers < 3 # ions are not present

    #     for icc in 1:data.params.numberOfCarriers # Integers
    #         if data.isContinuous[icc] == false
    #             data.chargeCarrierList[icc] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = icc)
    #         elseif data.isContinuous[icc] == true
    #             data.chargeCarrierList[icc] = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = icc)
    #         end
    #     end

    # else # ions are present
    #     ionic_carriers = data.enable_ionic_carriers.ionic_carriers

    #     for icc in 1:data.params.numberOfCarriers # Integers

    #         if data.isContinuous[icc] == false # discontinuous quantity
    #             if icc ∈ ionic_carriers # ionic quantity
    #                 data.chargeCarrierList[icc] = DiscontinuousQuantity(fvmsys, data.enable_ionic_carriers.regions, id = icc)
    #             else
    #                 data.chargeCarrierList[icc] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = icc)
    #             end
    #         end

    #     end

    # end

    for icc in 1:2# Integers
        if data.isContinuous[icc] == false
            data.chargeCarrierList[icc] = DiscontinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = icc)
        elseif data.isContinuous[icc] == true
            data.chargeCarrierList[icc] = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions, id = icc)
        end
    end

    if data.params.numberOfCarriers > 2
        data.chargeCarrierList[3] = InterfaceQuantity(fvmsys, 3, id = 3)
        data.chargeCarrierList[4] = InterfaceQuantity(fvmsys, 3, id = 4)
    end

    # put Auger and radiative on or off
    if data.bulk_recombination.bulk_recomb_Auger == false
        data.params.recombinationAuger .= 0.0
    end

    if data.bulk_recombination.bulk_recomb_radiative == false
        data.params.recombinationRadiative .= 0.0
    end

    # data.interfaceCarrierList = Array{VoronoiFVM.AbstractQuantity, 1}(undef, length(data.interfaceCarriers.index))
    # for ii = 1:length(data.interfaceCarriers.index)
    #     data.interfaceCarrierList[ii] = InterfaceQuantity(fvmsys, data.interfaceCarriers.boundary_region, id = data.interfaceCarriers.index[ii])
    # end

    data.index_psi        = ContinuousQuantity(fvmsys, 1:data.params.numberOfRegions)

    physics    = VoronoiFVM.Physics(data        = data,
                                    flux        = flux!,
                                    reaction    = reaction!,
                                    breaction   = breaction!,
                                    storage     = storage!,
                                    bstorage    = bstorage!,
                                    bflux       = bflux!
                                    )

    # add the defined physics to system
    physics!(fvmsys, physics)

    ctsys.fvmsys = fvmsys
    ctsys.data   = data

    return ctsys

end

###########################################################
###########################################################

function show_params(ctsys::System)

    params = ctsys.data.params
    for name in fieldnames(typeof(params))[1:end]
        @printf("%30s = ",name)
        println(getfield(params,name))
    end

end

function Base.show(io::IO, this::ParamsNodal)
    for name in fieldnames(typeof(this))[1:end]
        @printf("%30s = ",name)
        println(io,getfield(this,name))
    end
end

###########################################################
###########################################################

# Method which determines with input parameters which inner interface model was chosen by user.
function inner_interface_model(data::Data)

    countDiscontqF = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:data.params.numberOfBoundaryRegions

        if     data.boundary_type[ireg] == InterfaceModelDiscontqF

            countDiscontqF = countDiscontqF + 1

        elseif data.boundary_type[ireg] == InterfaceModelIonCharge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countDiscontqF > 0 # build discont_qF based system

        return InterfaceModelDiscontqF

    elseif countInterfaceCharge > 0 # build ion interface charge system

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in
        #     Example107 works perfectly fine since the choice of species number is already set.
        return InterfaceModelNone

    elseif countDiscontqF + countInterfaceCharge == 0 # build system without further interface conditions

        return InterfaceModelNone

    end

end


# Method which determines with input parameters which inner interface model was chosen by user.
function inner_interface_model(ctsys::System)


    countDiscontqF = 0::Int64; countInterfaceCharge = 0::Int64

    # detect which interface model the user chooses by counting
    for ireg in 1:ctsys.data.params.numberOfBoundaryRegions

        if     ctsys.data.boundary_type[ireg] ==  InterfaceModelDiscontqF

            countDiscontqF = countDiscontqF + 1

        elseif ctsys.data.boundary_type[ireg] == InterfaceModelIonCharge

            countInterfaceCharge = countInterfaceCharge + 1

        end

    end

     # build the system based on the input interface model
     if     countDiscontqF > 0 # build discont_qF based system

        return InterfaceModelDiscontqF

    elseif countInterfaceCharge > 0 # build ion interface charge system

        # DA: currently for this case, since InterfaceQuantites is not well tested we stick with the no inferface case, i.e.
        #     the known case how we build the system. Since we did not change anything, the code in
        #     Example107 works perfectly fine since the choice of species number is already set.
        return InterfaceModelNone

    elseif countDiscontqF + countInterfaceCharge == 0 # build system without further interface conditions

        return InterfaceModelNone

    end

end

###########################################################
###########################################################

"""
$(TYPEDSIGNATURES)

Master function which applies the voltage ``\\Delta u``at the
boundary ibreg for the chosen contact model.

"""

set_contact!(ctsys, ibreg, ;Δu) = __set_contact!(ctsys, ibreg, Δu, ctsys.data.boundary_type[ibreg])

# For schottky contacts
function __set_contact!(ctsys, ibreg, Δu, ::Type{SchottkyContact})

    ipsi = ctsys.data.index_psi

    # set Schottky barrier and applied voltage
    ctsys.fvmsys.boundary_values[ipsi,  ibreg] = (ctsys.data.params.SchottkyBarrier[ibreg]/q ) + Δu
    ctsys.fvmsys.boundary_factors[ipsi, ibreg] = VoronoiFVM.Dirichlet

end


function __set_contact!(ctsys, ibreg, Δu, ::Type{OhmicContact})

    interface_model       = ctsys.data.inner_interface_model

    set_ohmic_contact!(ctsys, ibreg, Δu, interface_model)

end

# DA: need to check, if the distinction here is necessary in future version
# (correlates with question of DiscontinuousQuantities)
function set_ohmic_contact!(ctsys, ibreg, Δu, ::Type{InterfaceModelNone})

    iphin = ctsys.data.bulk_recombination.iphin
    iphip = ctsys.data.bulk_recombination.iphip

    iphin = ctsys.data.chargeCarrierList[iphin]
    iphip = ctsys.data.chargeCarrierList[iphip]

    ctsys.fvmsys.boundary_factors[iphin, ibreg] = VoronoiFVM.Dirichlet
    ctsys.fvmsys.boundary_values[iphin, ibreg]  = Δu
    ctsys.fvmsys.boundary_factors[iphip, ibreg] = VoronoiFVM.Dirichlet
    ctsys.fvmsys.boundary_values[iphip, ibreg]  = Δu

end

function set_ohmic_contact!(ctsys, ibreg, contact_val, ::Type{InterfaceModelSurfaceReco})

    iphin = ctsys.data.bulk_recombination.iphin
    iphip = ctsys.data.bulk_recombination.iphip

    iphin = ctsys.data.chargeCarrierList[iphin]
    iphip = ctsys.data.chargeCarrierList[iphip]

    ctsys.fvmsys.boundary_factors[iphin, ibreg] = VoronoiFVM.Dirichlet
    ctsys.fvmsys.boundary_values[iphin, ibreg]  = contact_val
    ctsys.fvmsys.boundary_factors[iphip, ibreg] = VoronoiFVM.Dirichlet
    ctsys.fvmsys.boundary_values[iphip, ibreg]  = contact_val

end


function set_ohmic_contact!(ctsys, ibreg, contact_val, ::Type{InterfaceModelDiscontqF})

    iphin = ctsys.data.bulk_recombination.iphin
    iphip = ctsys.data.bulk_recombination.iphip

    iphin = ctsys.data.chargeCarrierList[iphin]
    iphip = ctsys.data.chargeCarrierList[iphip]

    if ibreg == 1

        for icc ∈ [iphin, iphip]
            ctsys.fvmsys.boundary_factors[icc.regionspec[ibreg], ibreg] = VoronoiFVM.Dirichlet
            ctsys.fvmsys.boundary_values[icc.regionspec[ibreg], ibreg] = contact_val
        end

    else

        for icc ∈ [iphin, iphip]
            ctsys.fvmsys.boundary_factors[icc.regionspec[ctsys.data.params.numberOfRegions], ibreg] = VoronoiFVM.Dirichlet
            ctsys.fvmsys.boundary_values[ icc.regionspec[ctsys.data.params.numberOfRegions], ibreg] = contact_val
        end

    end

end

###########################################################
###########################################################
# Wrappers for methods of VoronoiFVM

VoronoiFVM.enable_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_species!(ctsys.fvmsys, ispecies, regions)

VoronoiFVM.enable_boundary_species!(ctsys::System, ispecies, regions) = VoronoiFVM.enable_boundary_species!(ctsys.fvmsys, ispecies, regions)

VoronoiFVM.unknowns(ctsys::System) = VoronoiFVM.unknowns(ctsys.fvmsys)

VoronoiFVM.solve!(solution, initialGuess, ctsys, ;control=control, tstep=tstep) = VoronoiFVM.solve!(solution, initialGuess, ctsys.fvmsys, control=control, tstep=tstep)

###########################################################
###########################################################
"""
$(TYPEDSIGNATURES)

Functions which sets for given charge carrier at a given boundary
a given value.

"""

function equilibrium_solve!(ctsys::System; control = VoronoiFVM.NewtonControl(), nonlinear_steps = 20.0)

    ctsys.data.calculation_type    = InEquilibrium

    # initialize solution and starting vectors
    initialGuess                   = unknowns(ctsys)
    solution                       = unknowns(ctsys)
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

            VoronoiFVM.solve!(solution, initialGuess, ctsys, control = control, tstep=Inf)

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

"""
Gives the user a linear time mesh, this mesh is used for a linear I-V scan protocol.

"""
function set_time_mesh(scanrate, endVoltage, number_tsteps; type_protocol = LinearScanProtocol)

    tend  = endVoltage/scanrate

    return range(0, stop = tend, length = number_tsteps)
end


"""
Calculates current for time dependent problem.
"""
function get_current_val(ctsys, U, Uold, Δt) # DA: But caution, still need some small modification!

    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

    # left outer boundary = 1; right outer boundary = 2 (caution with order)
    tf     = testfunction(factory, [1], [2])
    I      = integrate(ctsys.fvmsys, tf, U, Uold, Δt)

    current = 0.0
    for icc in 1:ctsys.data.params.numberOfCarriers+1
        current = current + I[icc]
    end

    # DA: caution I[ipsi] not completly correct. In our examples, this does not effect something,
    # but we need derivative here.
    return current
end
###########################################################
###########################################################

"""
Calculates current for stationary problem.
"""
function get_current_val(ctsys, U) # DA: But caution, still need some small modification!

    factory = VoronoiFVM.TestFunctionFactory(ctsys.fvmsys)

    # left outer boundary = 1; right outer boundary = 2 (caution with order)
    tf     = testfunction(factory, [1], [2])
    I      = integrate(ctsys.fvmsys, tf, U)

    current = 0.0
    for icc in 1:ctsys.data.params.numberOfCarriers
        current = current + I[icc]
    end

    return current

end

###########################################################
###########################################################

"""

$(TYPEDSIGNATURES)

For given potentials, compute corresponding densities. This function is needed
for the method, plotting the densities.

"""
function compute_densities!(u, data, inode, region, icc, ipsi, in_region::Bool)

    params      = data.params
    paramsnodal = data.paramsnodal

    if in_region == false
        (params.bDensityOfStates[icc, region] + paramsnodal.densityOfStates[icc, inode] ) * data.F[icc](etaFunction(u, data, inode, region, icc, ipsi, in_region::Bool))
    elseif in_region == true
        (params.densityOfStates[icc, region] + paramsnodal.densityOfStates[icc, inode])* data.F[icc](etaFunction(u, data, inode, region, icc, ipsi, in_region::Bool))
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
    densities    = Array{Real,2}(undef, params.numberOfCarriers, size(sol, 2))

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

###########################################################
###########################################################

"""

$(SIGNATURES)

For given solution in vector form, compute corresponding vectorized band-edge energies and
Fermi level. [Caution: this was not tested for multidimensions.]
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
    regionsAllCells = push!(regionsAllCells, grid[CellRegions][end]) # enlarge region by final cell
    phi             = 0.0                                            # in equilibrium set to 0
    psi0_initial    = 0.5

    for index = 1:length(regionsAllCells) - 1

        ireg          = regionsAllCells[index]
        zVector       = params.chargeNumbers[iccVector]
        FVector       = data.F[iccVector]
        # all regions of nodes belonging to cell for given index
        regionsOfCell = regionsAllCells[grid[CellNodes][:,index]]

        # average following quantities if needed among all regions
        EVector = Float64[]; CVector = Float64[]; NVector = Float64[]

        for icc = 1:params.numberOfCarriers
            push!(EVector, sum(params.bandEdgeEnergy[icc, regionsOfCell])  / length(regionsOfCell) + paramsnodal.bandEdgeEnergy[icc,index])
            push!(CVector, sum(params.doping[icc, regionsOfCell])          / length(regionsOfCell) + paramsnodal.doping[index])
            push!(NVector, sum(params.densityOfStates[icc, regionsOfCell]) / length(regionsOfCell) + paramsnodal.densityOfStates[icc, index])
        end
        # rhs of Poisson's equation as anonymous function depending on psi0
        f = psi0 -> charge_density(psi0, phi, params.UT, EVector, zVector, CVector, NVector, FVector)

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
function charge_density(psi0, phi, UT, EVector, chargeNumbers, dopingVector, dosVector, FVector)
    # https://stackoverflow.com/questions/45667291/how-to-apply-one-argument-to-arrayfunction-1-element-wise-smartly-in-julia
    sum(-chargeNumbers .* dopingVector) + sum(chargeNumbers .* dosVector .* (etaFunction(psi0, phi, UT, EVector, chargeNumbers) .|> FVector))
end

"""

$(TYPEDSIGNATURES)

First try of debugger. Print the Jacobi matrix for a given node, i.e.
the number of the node in the grid and not the excact coordinate.
This is only done for the one dimensional case so far.
"""
function printJacobi(node, sys)
    ctdata = data(sys)
    numberOfNodes = ctdata.numberOfNodes
    if node == 1
        println(sys.matrix[1:3, 1:9])
    elseif node == numberOfNodes
        println(sys.matrix[3*numberOfNodes-2:3*numberOfNodes, 3*numberOfNodes-8:3*numberOfNodes])
    else
        println(sys.matrix[3*node-2:3*node, 3*node-5:3*node+3])
    end
end

"""

$(TYPEDSIGNATURES)

Compute the charge density for each region separately.
"""
function charge_density(ctsys, sol)
    integrate(ctsys.fvmsys,reaction!,sol)[ctsys.data.index_psi,:]
end