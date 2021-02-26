"""

Physical units in SI units

"""
# module Units
# using  Constants

# SI units
K                  = 1.0                      # Kelvin
J                  = 1.0                      # Joule
A                  = 1.0                      # Ampere
V                  = 1.0                      # Volt
m                  = 1.0                      # meter
s                  = 1.0                      # seconds
C                  = A * s                    # Coulomb
kg                 = 1.0                      # kilogramm

export K, J, A, V, m, s, C, kg

# scaled SI units
cm                 = 1.0e-2 * m
mm                 = 1.0e-3 * m
μm                 = 1.0e-6 * m 
nm                 = 1.0e-9 * m

ms                 = 1.0e-3 * s
μs                 = 1.0e-6 * s 
ns                 = 1.0e-9 * s
ps                 = 1.0e-12 *s

eV                 = q * V

export cm, mm, μm, nm, ms, μs, ns, ps, eV

# end
