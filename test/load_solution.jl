@testset "Load solution" begin

    filename = "notebook/NACA4415Re500.jld"
    u, t, sys = load_ns_solution(filename)

end
