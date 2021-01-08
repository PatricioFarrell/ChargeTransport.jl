using Documenter
using ChargeTransportInSolids

makedocs(sitename ="ChargeTransportInSolids.jl",
         modules  = [ChargeTransportInSolids],
         clean    = true,
         doctest  = false,
         authors  = "D. Abdel, P. Farrell, J. Fuhrmann",
         repo     = "https://github.com/PatricioFarrell/ChargeTransportInSolids.jl",
         pages    = [ "About the System" =>  "system.md",
                      "Plotting Functions" => "plot.md" ]
         )

# deploydocs(
#     repo = "github.com/PatricioFarrell/ChargeTransportInSolids.jl.git",
# )
