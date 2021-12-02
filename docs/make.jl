using Documenter
using ChargeTransportInSolids
using Literate


makedocs(sitename ="ChargeTransportInSolids.jl",
         modules  = [ChargeTransportInSolids],
         clean    = true,
         doctest  = false,
         authors  = "D. Abdel, P. Farrell, J. Fuhrmann",
         repo     = "https://github.com/PatricioFarrell/ChargeTransportInSolids.jl",
         pages    = Any[
            "ChargeTransportInSolids.jl -- A drift-diffusion solver" =>  "general.md",
            "Mathematical Description of the Problem" => "backgroundinfo.md",
            "Some Applications" => ["GeneralInformation.md",
                                    "examples/GaAs.md",
                                    "examples/PSC.md"],
            "Overview on Types and Constructors" => "allindex.md"
            ]
         )

deploydocs(
    repo = "github.com/PatricioFarrell/ChargeTransportInSolids.jl.git"
)
