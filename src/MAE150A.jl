module MAE150A

  using Reexport

  @reexport using ViscousFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq
  using JLD

  export save_ns_solution,load_ns_solution, compute_trajectory

  """
      save_ns_solution(filen::String,state,f,sys,body)

  Save a state vector `state` of the Navier-Stokes solution in a file with the provided
  filenam `filen`. Also save the associated vector of Lagrange forces `f`, the
  system solution metadata `sys`, and the body data `body`.
  """
  function save_ns_solution(filen,sys,body,state,f)

    save(filen,"Re",sys.Re,"U∞",sys.U∞,"Δt",sys.Δt,"grid",sys.grid,"body",body,"state",state,"f",f)
    return nothing

  end

  """
      load_ns_solution(filen::String)

  Given a file name `filen` in JLD format, load in the Navier-Stokes solution data
  stored in this file and set up various solution variables. An example:
  ```
  q, ω, ψ, f, sys, body = load_ns_solution("myfile.jld")
  ```
  In this example, `q` is the velocity field, `ω` the vorticity field,
  `ψ` the streamfunction, and `f` the vector of Lagrange forces on the body,
  `sys` is the NS system metadata and operators, and `body` the body data.
  """
  function load_ns_solution(filen)

    d = load(filen)

    body = d["body"]
    g = d["grid"]
    Δt = d["Δt"]
    U∞ = d["U∞"]
    Re = d["Re"]
    w = d["state"]
    f = d["f"]

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

    return q, ω, ψ, f, sys, body

  end

### TRAJECTORY CALCULATION ###

  function _vfcn!(dR,R,p,t,u,v)
    dR[1] = u(R[1],R[2])
    dR[2] = v(R[1],R[2])

   return dR
 end

 function compute_trajectory(ufield,vfield,X₀::Tuple,Tmax::Real,Δt::Real)

   u0 = [X₀[1],X₀[2]]
   tspan=(0.0,Tmax)

   vfcn!(dR,R,p,t) = _vfcn!(dR,R,p,t,ufield,vfield)

   Path = ODEProblem(vfcn!,u0,tspan)
   sol = solve(Path,ABM54(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)

   return sol
 end

end
