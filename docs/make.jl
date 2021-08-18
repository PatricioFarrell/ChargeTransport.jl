using Documenter
using ChargeTransportInSolids
using Literate


makedocs(sitename ="ChargeTransportInSolids.jl",
         modules  = [ChargeTransportInSolids],
         clean    = true,
         doctest  = false,
         authors  = "D. Abdel, P. Farrell, J. Fuhrmann",
         repo     = "https://github.com/PatricioFarrell/ChargeTransportInSolids.jl",
         pages    = [ "ChargeTransportInSolids.jl -- A drift-diffusion solver" =>  "general.md",
                      "Mathematical Description of the Problem" => "backgroundinfo.md",
                      "Numerical Strategy" =>  "system.md",
                      "Some Applications" => ["examples/GaAs.md", "examples/PSC.md"],
                      ##
                      "Plotting Functions" => "plot.md"
                      ]
         )

# deploydocs(
#     repo = "github.com/PatricioFarrell/ChargeTransportInSolids.jl.git",
# )
