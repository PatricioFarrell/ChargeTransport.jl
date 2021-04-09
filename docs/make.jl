using Documenter
using ChargeTransportInSolids
using Literate

# following function was just copied from VoronoiFVM
# Replace SOURCE_URL marker with github url of source
#
function replace_source_url(input,source_url)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out=IOBuffer()
    for line in lines_in
        println(lines_out,replace(line,"SOURCE_URL" => source_url))
    end
    return String(take!(lines_out))
end


example_jl_dir = joinpath(@__DIR__,"..","examples")
example_md_dir  = joinpath(@__DIR__,"src","examples")

for example_source in readdir(example_jl_dir) base,ext=splitext(example_source)
    if ext==".jl"
        source_url="https://github.com/PatricioFarrell/ChargeTransportInSolids.jl/tree/master/examples"*example_source
        preprocess(buffer)=replace_source_url(buffer,source_url)
        Literate.markdown(joinpath(@__DIR__,"..","examples",example_source),
                          example_md_dir,
                          documenter=false,
                          info=false,
                          preprocess=preprocess)
    end
end
generated_examples=vcat(["runexamples.md"],joinpath.("examples",readdir(example_md_dir)))

makedocs(sitename ="ChargeTransportInSolids.jl",
         modules  = [ChargeTransportInSolids],
         clean    = true,
         doctest  = false,
         authors  = "D. Abdel, P. Farrell, J. Fuhrmann",
         repo     = "https://github.com/PatricioFarrell/ChargeTransportInSolids.jl",
         pages    = [ "Introducing the Model" =>  "general.md",
                      "About the discrete System" =>  "system.md",
                      "Plotting Functions" => "plot.md",
                      "Examples" => generated_examples ]
         )

# deploydocs(
#     repo = "github.com/PatricioFarrell/ChargeTransportInSolids.jl.git",
# )
