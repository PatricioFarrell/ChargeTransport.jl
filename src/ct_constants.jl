
""" 
    Famous constants in SI units
"""

### Below is a first try with Unitful/PhysicalConstants. 
### Problem: Units need to be erased, see Quire #109 in ChargeTransportInSolids

# using Unitful
# using PhysicalConstants.CODATA2018

# kB                 = BoltzmannConstant                # JK^{-1}      --- Boltzmann constant
# Planck_constant    = PlanckConstant                   # Js           --- Planck constant
# mₑ                 = ElectronMass                     # kg           --- electron rest mass
# q                  = ElementaryCharge                 # C            --- elementary charge
# ε0                 = VacuumElectricPermittivity       # (A*s)/ (V*m) --- absolute dielectric permittivity of classical vacuum


const kB                 = 1.38064852e-23          # JK^{-1}      --- Boltzmann constant
const Planck_constant    = 6.62607015e-34          # Js           --- Planck constant
const mₑ                 = 9.1093837015e-31        # kg           --- electron rest mass
const q                  = 1.602176634e-19         # C            --- elementary charge
const ε0                 = 8.8541878176e-12        # (A*s)/ (V*m) --- absolute dielectric permittivity of classical vacuum


# Set this as temporary values from pdelib, figure out the recent SI stuff (Unitful.jl ?
function set_pdelib_constants()
    global kB = 1.3806503e-23
    global q  = 1.602176462e-19
    global ε0 = 8.85418781762039e-12
    global eV = q * V
end

function set_TeSCA_constants()
    global kB = 1.380662e-23
    global q  = 1.6021e-19
    global ε0 = 8.85419e-12
    global eV = q * V
end
