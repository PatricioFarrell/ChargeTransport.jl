

# Default parameters of Ionmonger (https://github.com/PerovskiteSCModelling/IonMonger)
# representing TiO2 | MAPI | spiro-OMeTAD

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
h_ndoping          = 100.0 * nm
h_intrinsic        = 400.0 * nm
h_pdoping          = 200.0 * nm
h_total            = h_ndoping + h_intrinsic + h_pdoping
heightLayers       = [h_ndoping,
                            h_ndoping + h_intrinsic,
                            h_ndoping + h_intrinsic + h_pdoping]


########## physical values ##########

## charge numbers
zn                 = -1
zp                 = 1
za                 = 1

## temperature
T                  = 298.0                          *  K

## band edge energies
En                 = [-4.0, -3.7, -3.4]            .*  eV
Ep                 = [-5.8, -5.4, -5.1]            .*  eV
Ea                 = [0.0, -4.45,  0.0]            .*  eV
Ea_i               = Ea[regionIntrinsic]

## effective densities of density of states
Nn                 = [5.0e25, 8.1e24, 5.0e25]      ./ (m^3)
Np                 = [5.0e25, 5.8e24, 5.0e25]      ./ (m^3)
Na                 = [0.0,    1.0e27, 0.0]         ./ (m^3)
Na_i               = Na[regionIntrinsic]

## mobilities
μn                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
μp                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
μa                 = [0.0, 3.93e-16, 0.0]          .* (m^2) / (V * s)

## relative dielectric permittivity
ε                  = [10.0, 24.1, 3.0]             .* 1.0

## radiative recombination
r0                 = [6.8e-17, 3.6e-18, 6.3e-17]   .* m^3 / s

## life times and trap densities
τn                 = [1.0e100, 3.0e-9, 1.0e100]    .* s
τp                 = [1.0e100, 3.0e-7, 1.0e100]    .* s

## SRH trap energies
EI                 = [-5.0, -4.55, -4.1]           .* eV

## generation
incidentPhotonFlux = [0.0, 9.4e20, 0.0]            ./ (m^2 * s)
absorption         = [0.0, 1.3e7, 0.0]             ./ m
generationPeak     = h_ndoping

generation_uniform = [0.0, 2.5e27, 0.0]           ./ (m^3 * s)

## doping
Cn                 = 1.0e24                         / (m^3)
Cp                 = 1.0e24                         / (m^3)
Ca                 = 1.6e25                         / (m^3)
