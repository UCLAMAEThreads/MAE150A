using MAE150A
using ViscousFlow
using Test
using Literate

#include("trajectories.jl")

outputdir = "../outputtest"
for (root, dirs, files) in walkdir("literate")
    if splitpath(root)[end] == "assets"
        for file in files
            cp(joinpath(root, file),joinpath(outputdir,file),force=true)
        end
    end
end

for (root, dirs, files) in walkdir("literate")
    for file in files
        endswith(file,".jl") && Literate.notebook(joinpath(root, file),outputdir=outputdir)
    end
end
