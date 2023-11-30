

# Default parameters of Ionmonger (https://github.com/PerovskiteSCModelling/IonMonger)
# representing TiO2 | MAPI | spiro-OMeTAD

#####################################################################
############################ parameters ############################

########## charge carriers ##########

const iphin              = 1 # electron quasi Fermi potential
const iphip              = 2 # hole quasi Fermi potential
const iphia              = 3
const numberOfCarriers   = 3 # electrons, holes and anion vacancies

########## device geometry ##########

# region numbers
const regionDonor        = 1
const regionIntrinsic    = 2
const regionAcceptor     = 3
const regions            = [regionDonor, regionIntrinsic, regionAcceptor]
const numberOfRegions    = length(regions)

# boundary region numbers
const bregionDonor       = 1
const bregionAcceptor    = 2
const bregionJ1          = 3
const bregionJ2          = 4

## length domains
const h_ndoping          = 100.0 * nm
const h_intrinsic        = 400.0 * nm
const h_pdoping          = 200.0 * nm
const h_total            = h_ndoping + h_intrinsic + h_pdoping
const heightLayers       = [h_ndoping,
                            h_ndoping + h_intrinsic,
                            h_ndoping + h_intrinsic + h_pdoping]


########## physical values ##########

## charge numbers
const zn                 = -1
const zp                 = 1
const za                 = 1

## temperature
const T                  = 298.0                          *  K

## band edge energies
const En                 = [-4.0, -3.7, -3.4]            .*  eV
const Ep                 = [-5.8, -5.4, -5.1]            .*  eV
const Ea                 = [0.0, -4.45,  0.0]            .*  eV
const Ea_i               = Ea[regionIntrinsic]

## effective densities of density of states
const Nn                 = [5.0e25, 8.1e24, 5.0e25]      ./ (m^3)
const Np                 = [5.0e25, 5.8e24, 5.0e25]      ./ (m^3)
const Na                 = [0.0,    1.0e27, 0.0]         ./ (m^3)
const Na_i               = Na[regionIntrinsic]

## mobilities
const μn                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
const μp                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
const μa                 = [0.0, 3.93e-16, 0.0]          .* (m^2) / (V * s)

## relative dielectric permittivity
const ε                  = [10.0, 24.1, 3.0]             .* 1.0

## radiative recombination
const r0                 = [6.8e-17, 3.6e-18, 6.3e-17]   .* m^3 / s

## life times and trap densities
const τn                 = [1.0e100, 3.0e-9, 1.0e100]    .* s
const τp                 = [1.0e100, 3.0e-7, 1.0e100]    .* s

## SRH trap energies
const EI                 = [-5.0, -4.55, -4.1]           .* eV

## generation
const incidentPhotonFlux = [0.0, 9.4e20, 0.0]            ./ (m^2 * s)
const absorption         = [0.0, 1.3e7, 0.0]             ./ m
const generationPeak     = h_ndoping

const generation_uniform = [0.0, 2.64e27, 0.0]           ./ (m^3 * s)

## doping
const Cn                 = 1.0e24                         / (m^3)
const Cp                 = 1.0e24                         / (m^3)
const Ca                 = 1.6e25                         / (m^3)
