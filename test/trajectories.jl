using Plots
using ViscousFlow

@testset "Trajectories" begin

    filename = "../notebook/NACA4415Re500.jld"
    u, t, sys = load_ns_solution(filename)

    vel = ViscousFlow.velocity(u,sys,t);


    pts = [ [-1,0.25], [-1,0], [-1,-0.025], [-1,-0.05], [-1,-0.25],[-1,-0.5],  [-1,-0.75] ]

    Tmax = 10.0

    traj_array = compute_trajectory(vel,sys,pts,Tmax)

    trajectories(traj_array,sys)

end
