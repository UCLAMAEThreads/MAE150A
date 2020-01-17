module MAE150A

  using Reexport

  @reexport using ViscousFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq
  @reexport using LaTeXStrings
  using Interpolations
  using JLD
  #using PyPlot

  export initialize_environment,initialize_ns_solver,
        save_ns_solution,load_ns_solution, get_flowfield,
        compute_trajectory, field_along_trajectory,
        convective_acceleration, ddt

  function initialize_environment()

    # Set the back end for Plots
    #pyplot()
    rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

    # Ensure that LaTeX stuff is handled
    rcParams["text.usetex"] = true

    return nothing

  end

  """
      initialize_ns_solver(Re::Real,U∞::Tuple,Δx::Real,xlim::Tuple,ylim::Tuple,body::Body[,Δt = Nothing])

  Initialize the Navier-Stokes solver for a problem, with Reynolds number `Re`,
  free stream velocity `U∞` on a grid with cell size `Δx` and extent `xlim`
  and `ylim`, and body `body`. The time step size is set automatically, unless it
  is passed in as an optional argument.

  An example usage is

  ```
  solver,sys,state,f = initialize_ns_solver(Re,U∞,Δx,xlim,ylim,body)
  ```

  The output `solver` is the integrator to use for advancing the system, `sys`
  contains operators and the NS system's metadata, `state` is an example instance of the
  state vector, and `f` and example instance of the constraint force vector.
  """
  function initialize_ns_solver(Re,U∞,Δx,xlim,ylim,body::Body;Δt = Nothing)
      if Δt == Nothing
          # Use CFL criteria to establish time step size
          Δt = min(0.5*Δx,0.5*Δx^2*Re)
      end
      X = VectorData(body.x,body.y)

      sys = NavierStokes(Re,Δx,xlim,ylim,Δt,U∞ = U∞, X̃ = X, isstore = true)
      state = Nodes(Dual,size(sys))
      f = VectorData(X)

      plan_intfact(t,u) = Systems.plan_intfact(t,u,sys)
      plan_constraints(u,t) = TimeMarching.plan_constraints(u,t,sys)
      r₁(u,t) = TimeMarching.r₁(u,t,sys)
      r₂(u,t) = TimeMarching.r₂(u,t,sys)

      return IFHERK(state,f,sys.Δt,plan_intfact,plan_constraints,(r₁,r₂),rk=TimeMarching.RK31,isstored=true), sys, state, f
  end


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
      get_flowfield(state,sys::NavierStokes)

  Get the other flow field quantities: velocity, vorticity, and streamfunction
  from the given state `state`, for Navier-Stokes system `sys`. Usage:

  ```
  u, ω, ψ = get_flowfield(state,sys)
  ```
  """
  function get_flowfield(w::Nodes{Dual},sys)
    xg, yg = coordinates(w,sys.grid)
    ω = vorticity(w,sys)
    q = velocity(w,sys)

    q.u .+= sys.U∞[1]
    q.v .+= sys.U∞[2]
    ψ = streamfunction(w,sys)
    ψ .+= sys.U∞[1]*transpose(yg)

    return q, ω, ψ
  end

  """
      load_ns_solution(filen::String)

  Given a JLD file with name `filen`, load in the Navier-Stokes solution data
  stored in this file and set up various solution variables. An example:

  ```
  state, f, sys, body = load_ns_solution("myfile.jld")
  ```

  In this example, `state` is the flow state vector, `f` the vector of Lagrange
  forces on the body, `sys` is the NS system metadata and operators, and `body` the body data.
  """
  function load_ns_solution(filen)

    d = load(filen)

    body = d["body"]
    g = d["grid"]
    Δt = d["Δt"]
    U∞ = d["U∞"]
    Re = d["Re"]
    state = d["state"]
    f = d["f"]

    X = VectorData(body.x,body.y)
    sys = NavierStokes(Re,cellsize(g),limits(g,1),limits(g,2),Δt,
                        U∞ = U∞, X̃ = X, isstore = true)

    #q, ω, ψ = get_flowfield(w,sys)

    #return q, ω, ψ, f, sys, body
    return state, f, sys, body
  end

  # TRAJECTORY CALCULATION #

  """
     compute_trajectory(u,v,X₀::Vector/Vector{Vector},Tmax,Δt)

  Calculate the trajectory of a particle with initial location(s) `X₀`. More than
  one particle can be tracked by providing a vector of locations (a vector of 2-
  dimensional vectors) in X₀. The arguments
  `u` and `v` are interpolated velocity field components, `Tmax` is the final
  integration time, and `Δt` is the time step size. The output is the solution
  structure for the `OrdinaryDiffEq` package (or, for multiple particles, a vector
  of such solution structures).
  """
  function compute_trajectory(ufield::AbstractInterpolation{T,2},
                              vfield::AbstractInterpolation{T,2},
                              X₀::Vector{S},Tmax::Real,Δt::Real) where {T,S<:Real}

    u0 = X₀
    tspan=(0.0,Tmax)

    vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

    Path = ODEProblem(vfcn!,u0,tspan)
    sol = solve(Path,ABM54(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)

    return sol
  end

  function compute_trajectory(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},pts::Vector{Vector{S}},Tmax,Δt) where {T,S<:Real}

    sol_array = ODESolution[]
    for X₀ in pts
      sol = compute_trajectory(ufield,vfield,X₀,Tmax,Δt)
      push!(sol_array,sol)
    end
    return sol_array

  end

  """
     compute_trajectory(vel::Edges,sys,X₀::Vector/Vector{Vector},Tmax,Δt)

  Calculate the trajectory of a particle with initial location(s) `X₀`. The argument
  `vel` is edge-type grid data, `sys` is a Navier-Stokes type system, `Tmax` is the final
  integration time, and `Δt` is the time step size. The output is the solution
  structure for the `OrdinaryDiffEq` package.
  """
  compute_trajectory(u::Edges, sys::NavierStokes, X₀,Tmax,Δt) =
      compute_trajectory(interpolatable_field(u,sys.grid)...,X₀,Tmax,Δt)


  function _vfcn_autonomous!(dR,R,p,t,u,v)
    dR[1] = u(R[1],R[2])
    dR[2] = v(R[1],R[2])

   return dR
 end

 """
    field_along_trajectory(f::GridData,sys::NavierStokes,traj::ODESolution)

 Evaluate field `f` (given as grid data) along the trajectory specified by `traj`.
 The output is the history of `f` along this trajectory. If `f` is a vector field,
 then the component histories are output as a tuple.
 """
 function field_along_trajectory(v::VectorGridData,sys::NavierStokes,traj::ODESolution)
   vfield_x, vfield_y = interpolatable_field(v,sys.grid)

   vx_traj = eltype(v)[]
   vy_traj = eltype(v)[]
   for x in traj.u
     push!(vx_traj,vfield_x(x...))
     push!(vy_traj,vfield_y(x...))
   end

   return vx_traj, vy_traj
 end

 function field_along_trajectory(s::ScalarGridData,sys::NavierStokes,traj::ODESolution)
   sfield = interpolatable_field(s,sys.grid)

   s_traj = eltype(sfield)[]
   for x in traj.u
     push!(s_traj,sfield(x...))
   end

   return s_traj
 end

 # Convective acceleration
 """
    convective_acceleration(u,sys)

 Given grid velocity data `u` for Navier-Stokes system `sys`, calculation
 the convective acceleration field on the grid u.grad(u).
 """
 function convective_acceleration(u::VectorGridData,sys::NavierStokes)
   ugradu = zero(u)
   convective_derivative!(ugradu,u)
   ugradu ./= cellsize(sys.grid)

   return ugradu

 end

 """
    ddt(u::AbstractVector,Δt[,mydiff=:backward_diff])

 Calculate the time derivative of vector data `u`, with time step size `Δt`.
 The default method is backward differencing, but this can be changed to
 `:forward_diff` or `:central_diff`.
 """
 function ddt(u::AbstractVector{T},Δt::Real;mydiff::Symbol=:forward_diff) where {T}
    return eval(mydiff)(u)/Δt
 end

 # Some basic differencing routines
 function backward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end] .= u[2:end] .- u[1:end-1]
     u[1] = u[2]
     return du
 end
 function forward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[1:end-1] .= u[2:end] .- u[1:end-1]
     u[end] = u[end-1]
     return du
 end
 function central_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end-1] .= 0.5*u[3:end] .- 0.5*u[1:end-2]
     u[1] = u[2]
     u[end] = u[end-1]
     return du
 end



end
