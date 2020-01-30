module MAE150A

  using Reexport

  @reexport using ViscousFlow
  @reexport using PotentialFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq
  @reexport using LaTeXStrings

  using Interpolations
  using JLD
  #using PyPlot

  export initialize_environment,initialize_ns_solver,
        save_ns_solution,load_ns_solution, get_flowfield,
        compute_trajectory, field_along_trajectory,
        convective_acceleration, mag, ddt, pressure,
        OseenVortex,
        complexgrid, vortex_patch

  function initialize_environment()

    # Set the back end for Plots
    #pyplot()
    rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

    # Ensure that LaTeX stuff is handled
    rcParams["mathtext.fontset"] = "cm"
    #=
    This does not always work well...
    if typeof(Plots.PyPlot.matplotlib.checkdep_dvipng()) != Nothing
      # only use matplotlib tex if dvipng is present
      rcParams["text.usetex"] = true
    end
    =#

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
      get_flowfield(state,f,sys::NavierStokes)

  Get the other flow field quantities: velocity, vorticity, streamfunction,
  and pressure coefficient (Cp), from the given state `state` and associated
  constraint force vector `f`, for Navier-Stokes system `sys`. Usage:

  ```
  u, ω, ψ, Cp = get_flowfield(state,f,sys)
  ```
  """
  function get_flowfield(w::Nodes{Dual},f::VectorData,sys)
    xg, yg = coordinates(w,sys.grid)
    ω = vorticity(w,sys)
    q = velocity(w,sys)

    q.u .+= sys.U∞[1]
    q.v .+= sys.U∞[2]
    ψ = ViscousFlow.streamfunction(w,sys)
    ψ .+= sys.U∞[1]*transpose(yg)

    Cp = pressure(w,f,sys)

    return q, ω, ψ, Cp
  end

  function pressure(w::Nodes{Dual},f::VectorData,sys)

    u = velocity(w,sys)
    u.u .+= sys.U∞[1]
    u.v .+= sys.U∞[2]

    u_dual = Nodes(Dual,u)
    ucrossw = Edges(Primal,u)

    grid_interpolate!(ucrossw.u,grid_interpolate!(u_dual, u.v) ∘ w)
    grid_interpolate!(ucrossw.v,grid_interpolate!(u_dual,-u.u) ∘ w)
    rhs = divergence(-cellsize(sys)*(sys.Hmat*f) + ucrossw)

    umag = mag(u)

    Lc = plan_laplacian(rhs,with_inverse=true)

    fact = 2/(sys.U∞[1]^2+sys.U∞[2]^2)
    Cp = fact*(Lc\rhs - 0.5*(umag∘umag))

    return Cp

  end

  function pressure(w::Nodes{Dual},g::PhysicalGrid; U∞::Tuple=(0,0))

    L = plan_laplacian(w,with_inverse=true)

    u = -curl(L\w)
    u.u .+= U∞[1]
    u.v .+= U∞[2]

    u_dual = Nodes(Dual,u)
    ucrossw = Edges(Primal,u)

    grid_interpolate!(ucrossw.u,grid_interpolate!(u_dual, u.v) ∘ w)
    grid_interpolate!(ucrossw.v,grid_interpolate!(u_dual,-u.u) ∘ w)
    rhs = divergence(ucrossw)

    umag = mag(u)

    Lc = plan_laplacian(rhs,with_inverse=true)

    return Lc\rhs - 0.5*(umag∘umag)

  end

# Vortex construction

  """
      OseenVortex(x0,y0,Γ,σ)

  Construct a Gaussian-shaped vorticity distribution at location `x0`, `y0`, with
    strength `Γ` and radius `σ`. This is used as in the following example:
  ```
  w = OseenVortex(0,1,1,0.1)
  ```
  constructs a vortex at (0,1) with strength 1 and radius 0.1. Then we can evaluate
  the vorticity field of this vortex at some location (x,y) as follows:
  ```
  w(x,y)
  ```
  """
  struct OseenVortex
    x0 :: Real
    y0 :: Real
    Γ :: Real
    σ :: Real
  end

  (v::OseenVortex)(x,y) = v.Γ/(π*v.σ^2)*exp(-((x-v.x0)^2+(y-v.y0)^2)/v.σ^2)

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

 # Vector data magnitude (on cell centers)
 """
    mag(u::Edges{Primal/Dual}) -> Nodes{Primal/Dual}

 Calculate the magnitude of vector grid data `u`, placing the result on
 the cell centers.
 """
 function mag(u::Edges{C}) where {C <: ViscousFlow.Fields.CellType}

   usq = u∘u
   usq_nodes = Nodes(C,u)
   umag = Nodes(C,u)

   grid_interpolate!(usq_nodes,usq.u)
   umag .= usq_nodes

   grid_interpolate!(usq_nodes,usq.v)
   umag .+= usq_nodes

   @. umag = sqrt(umag)

   return umag
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

## Potential flow routines

function complexgrid(x::AbstractVector,y::AbstractVector)
    z = zeros(ComplexF64,length(x),length(y))
    @. z = x + im*y'
    return z
end

"""
    vortex_patch(xcent,ycent,strength,radius,nring) -> Vector{Vortex.Point}

Create a list of point vortices in the form of a vortex patch, a set of `nring` concentric
rings centered at `(xcent,ycent)` with radius `radius`. The overall circulation of the
patch is defined by `strength`.
"""
function vortex_patch(xcent,ycent,strength,radius,nring)
    Δr = radius/(nring-1)
    zcent = xcent+im*ycent

    r = 0.0

    zvort = ComplexF64[]
    cnt = 0
    for i = 1:nring
        θ = 0.0
        nv = max(1,8*(i-1))
        r = (i-1)*Δr
        for j = 1:nv
            push!(zvort,zcent + r*exp(im*θ))
            cnt += 1
            θ += 2π/nv
        end
    end

    return Vortex.Point.(zvort,strength/cnt)

end


end
