using Documenter
using ChargeTransport
using Literate


makedocs(sitename ="ChargeTransport.jl",
         modules  = [ChargeTransport],
         clean    = true,
         doctest  = false,
         authors  = "D. Abdel, P. Farrell, J. Fuhrmann",
         repo     = "https://github.com/PatricioFarrell/ChargeTransport.jl",
         pages    = Any[
            "ChargeTransport.jl -- A drift-diffusion solver" =>  "general.md",
            "Mathematical Description of the Problem" => "backgroundinfo.md",
            "Some Applications" => ["GeneralInformation.md",
                                    "examples/GaAs.md",
                                    "examples/PSC.md"],
            "Overview on Types and Constructors" => "allindex.md"
            ]
         )

deploydocs(
    repo = "github.com/PatricioFarrell/ChargeTransport.jl.git"
)
