module MAE150A

  using Reexport

  @reexport using ViscousFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq
  using JLD


  function load_ns_solution(filen)

    d = load(filen)

    body = d["body"]
    g = d["grid"]
    Δt = d["Δt"]
    U∞ = d["U∞"]
    Re = d["Re"]
    w = d["u"]

    X = VectorData(body.x,body.y);
    sys = NavierStokes(Re,cellsize(g),limits(g,1),limits(g,2),Δt,
                        U∞ = U∞, X̃ = X, isstore = true)

    xg, yg = coordinates(w,sys.grid)
    ω = vorticity(w,sys)
    q = velocity(w,sys)

    q.u .+= sys.U∞[1]
    q.v .+= sys.U∞[2]
    ψ = streamfunction(w,sys) #.+ sys.U∞[1]*yg'
    ψ .+= sys.U∞[1]*yg';

    return sys, body, q, ω, ψ

  end

  function _vfcn!(dR,R,p,t,u,v)
    dR[1] = u(R[1],R[2])
    dR[2] = v(R[1],R[2])

   return dR
 end

 vfcn!(dR,R,p,t) = _vfcn!(dR,R,p,t,ufield,vfield)

 function compute_trajectory(ufield,vfield,X₀::Tuple,Tmax::Real,Δt::Real)

   u0 = [X₀[1],X₀[2]]
   tspan=(0.0,Tmax)

   vfcn!(dR,R,p,t) = _vfcn!(dR,R,p,t,ufield,vfield)

   Path = ODEProblem(vfcn!,u0,tspan)
   sol = solve(Path,ABM54(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)

   return sol
 end

end
