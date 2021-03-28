using MAE150A
using ViscousFlow
using Test
using Literate

#include("trajectories.jl")

outputdir = "../notebook"
litdir = "./literate"

for (root, dirs, files) in walkdir(litdir)
    if splitpath(root)[end] == "assets"
        for file in files
            cp(joinpath(root, file),joinpath(outputdir,file),force=true)
        end
    end
end

function replace_includes(str)

    included = ["header.jl"]

    # Here the path loads the files from their proper directory,
    # which may not be the directory of the `examples.jl` file!
    path = litdir

    for ex in included
        content = read(joinpath(litdir,ex), String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

for (root, dirs, files) in walkdir(litdir)
    for file in files
        #endswith(file,".jl") && startswith(file,"4") && Literate.notebook(joinpath(root, file),outputdir,preprocess = replace_includes)
        endswith(file,".jl") && Literate.notebook(joinpath(root, file),outputdir,preprocess = replace_includes)
    end
end
