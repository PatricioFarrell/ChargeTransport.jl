

# Parameters from Driftfusion:
# (https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv)

# representing PCBM | MAPI | Pedot:PSS

#####################################################################
############################ parameters ############################

########## charge carriers ##########

iphin              = 1 # electron quasi Fermi potential
iphip              = 2 # hole quasi Fermi potential
iphia              = 3
numberOfCarriers   = 3 # electrons, holes and anion vacancies

########## device geometry ##########

# region numbers
regionDonor        = 1
regionIntrinsic    = 2
regionAcceptor     = 3
regions            = [regionDonor, regionIntrinsic, regionAcceptor]
numberOfRegions    = length(regions)

# boundary region numbers
bregionDonor       = 1
bregionAcceptor    = 2
bregionJ1          = 3
bregionJ2          = 4

## length domains
h_ndoping          = 60.0  * nm
h_intrinsic        = 300.0 * nm
h_pdoping          = 50.0  * nm
h_total            = h_ndoping + h_intrinsic + h_pdoping
const heightLayers       = [h_ndoping,
                            h_ndoping + h_intrinsic,
                            h_ndoping + h_intrinsic + h_pdoping]
########## physical values ##########

## charge numbers
const zn                 = -1
const zp                 = 1
const za                 = 1

## temperature
const T                  = 300.0                           *  K

## band edge energies
const En                 = [-3.8, -3.8,  -3.0]            .*  eV
const Ep                 = [-6.2, -5.4,  -5.1]            .*  eV
const Ea                 = [0.0, -4.66,  0.0]             .*  eV

## effective densities of density of states
const Nn                 = [1.0e25, 1.0e25,  1.0e26]      ./ (m^3)
const Np                 = [1.0e25, 1.0e25,  1.0e26]      ./ (m^3)
const Na                 = [0.0,    1.0e26, 0.0]          ./ (m^3)


## mobilities
const μn                 = [1.0e-7, 2.0e-3, 1.0e-5]       .* (m^2) / (V * s)
const μp                 = [1.0e-7, 2.0e-3, 1.0e-5]       .* (m^2) / (V * s)
const μa                 = [0.0, 1.0e-14,  0.0]           .* (m^2) / (V * s)


## relative dielectric permittivity
const ε                  = [3.0, 23.0, 4.0]               .* 1.0

## radiative recombination
const r0                 = [6.8e-17, 3.6e-18, 6.3e-17]    .* cm^3 / s

## life times and trap densities
const τn                 = [1.0e-6, 1.0e-7, 1.0e-6]       .* s
const τp                 = [1.0e-6, 1.0e-7, 1.0e-6]       .* s

## SRH trap energies
const EI                 = [-5.0, -4.60, -4.05]          .* eV

## generation
const incidentPhotonFlux = [0.0, 8.0e20, 0.0]             ./ (m^2 * s)
const absorption         = [0.0, 1.3e7, 0.0]              ./ m
const generationPeak     = h_ndoping

const generation_uniform = [0.0, 2.64e27, 0.0]            ./ (m^3 * s)

## doping
const Cn                 = 2.09e24                         / (m^3)
const Cp                 = 2.09e24                         / (m^3)
const Ca                 = 1.0e24                          / (m^3)

const UT                 = kB * T / q
