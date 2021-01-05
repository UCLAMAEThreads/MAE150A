module MAE150A

  using Weave, Pkg, InteractiveUtils, IJulia

  using Reexport

  @reexport using ViscousFlow
  @reexport using PotentialFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq
  @reexport using LaTeXStrings

  using Interpolations
  using JLD
  using PyCall
  using Dierckx
  using Roots
  #using PyPlot

  export initialize_environment,initialize_ns_solver,
        save_ns_solution,load_ns_solution, get_flowfield,
        compute_trajectory, compute_trajectories, field_along_trajectory,
        convective_acceleration, mag, ddt, pressure,
        OseenVortex,
        complexgrid, vortex_patch,
        add_arrow!,add_arrows!,
        falknerskan


  repo_directory = joinpath(@__DIR__,"..")
  cssfile = joinpath(@__DIR__, "..", "templates", "skeleton_css.css")
  latexfile = joinpath(@__DIR__, "..", "templates", "julia_tex.tpl")

  function weave_file(folder,file,build_list=(:script,:html,:pdf,:github,:notebook); kwargs...)
    tmp = joinpath(repo_directory,"tutorials",folder,file)
    Pkg.activate(dirname(tmp))
    Pkg.instantiate()
    args = Dict{Symbol,String}(:folder=>folder,:file=>file)
    if :script ∈ build_list
      println("Building Script")
      dir = joinpath(repo_directory,"script",folder)
      isdir(dir) || mkdir(dir)
      args[:doctype] = "script"
      tangle(tmp;out_path=dir)
    end
    if :html ∈ build_list
      println("Building HTML")
      dir = joinpath(repo_directory,"html",folder)
      isdir(dir) || mkdir(dir)
      args[:doctype] = "html"
      weave(tmp,doctype = "md2html",out_path=dir,args=args; fig_ext=".svg", css=cssfile, kwargs...)
    end
    if :pdf ∈ build_list
      println("Building PDF")
      dir = joinpath(repo_directory,"pdf",folder)
      isdir(dir) || mkdir(dir)
      args[:doctype] = "pdf"
      try
        weave(tmp,doctype="md2pdf",out_path=dir,args=args; template=latexfile, kwargs...)
      catch ex
        @warn "PDF generation failed" exception=(ex, catch_backtrace())
      end
    end
    if :github ∈ build_list
      println("Building Github Markdown")
      dir = joinpath(repo_directory,"markdown",folder)
      isdir(dir) || mkdir(dir)
      args[:doctype] = "github"
      weave(tmp,doctype = "github",out_path=dir,args=args; kwargs...)
    end
    if :notebook ∈ build_list
      println("Building Notebook")
      dir = joinpath(repo_directory,"notebook",folder)
      isdir(dir) || mkdir(dir)
      args[:doctype] = "notebook"
      Weave.convert_doc(tmp,joinpath(dir,file[1:end-4]*".ipynb"))
    end
  end

  function weave_all()
    for folder in readdir(joinpath(repo_directory,"tutorials"))
      folder == "test.jmd" && continue
      weave_folder(folder)
    end
  end

  function weave_folder(folder)
    for file in readdir(joinpath(repo_directory,"tutorials",folder))
      println("Building $(joinpath(folder,file)))")
      try
        weave_file(folder,file)
      catch
      end
    end
  end

  function tutorial_footer(folder=nothing, file=nothing; remove_homedir=true)
      display("text/markdown", """
      ## Appendix
       This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
       For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.
      """)
      if folder !== nothing && file !== nothing
          display("text/markdown", """
          To locally run this tutorial, do the following commands:
          ```
          using SciMLTutorials
          SciMLTutorials.weave_file("$folder","$file")
          ```
          """)
      end
      display("text/markdown", "Computer Information:")
      vinfo = sprint(InteractiveUtils.versioninfo)
      display("text/markdown",  """
      ```
      $(vinfo)
      ```
      """)

      ctx = Pkg.API.Context()
      pkgs = Pkg.Display.status(Pkg.API.Context(), use_as_api=true);
      projfile = ctx.env.project_file
      remove_homedir && (projfile = replace(projfile, homedir() => "~"))

      display("text/markdown","""
      Package Information:
      """)

      md = ""
      md *= "```\nStatus `$(projfile)`\n"

      for pkg in pkgs
          if !isnothing(pkg.old) && pkg.old.ver !== nothing
            md *= "[$(string(pkg.uuid))] $(string(pkg.name)) $(string(pkg.old.ver))\n"
          else
            md *= "[$(string(pkg.uuid))] $(string(pkg.name))\n"
          end
      end
      md *= "```"
      display("text/markdown", md)
  end

  function open_notebooks()
    Base.eval(Main, Meta.parse("import IJulia"))
    path = joinpath(repo_directory,"notebook")
    IJulia.notebook(;dir=path)
  end


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

    default(markerstrokealpha = 0, legend = false,
        dpi = 100, size = (400, 300), grid = false)

    return nothing

  end

### Plotting arrows on streamlines

"""
    get_segment_coords(p::Plot,level::Integer)

Get the x,y coordinates of contour level `level` in plot `p`.
"""
function get_segment_coords(ps::Plots.Plot{Plots.PyPlotBackend},level::Integer)
    contourset = get_series_handle(ps)

    if level > length(contourset."collections")
      error("This contour level is not on the list.")
    end
    lc = get(contourset."collections",level-1)
    seg_array_list = pycall(lc."get_segments",Array)

    return seg_array_list
end

function get_series_handle(ps::Plots.Plot{Plots.PyPlotBackend})
    sl = ps.series_list[1]
    contourset = sl.plotattributes[:serieshandle][1]
    return contourset
end

function arrowhead_coords(x,y,u,v,scale;open_angle::Float64=π/3)
    R = [cos(0.5*open_angle) sin(0.5*open_angle);
          -sin(0.5*open_angle) cos(0.5*open_angle)]

    uvec = [u,v]/sqrt(u^2+v^2)

    xhead = [x,y]
    x1 = xhead .- scale*(R*uvec)
    x2 = xhead .- scale*(R'*uvec)

    xdata = [x1[1],x,x2[1]]
    ydata = [x1[2],y,x2[2]]
    return xdata, ydata
end

function arrowhead!(ps::Plots.Plot{Plots.PyPlotBackend},x,y,u,v,scale;open_angle::Float64=π/3)

    xdata, ydata = arrowhead_coords(x,y,u,v,scale,open_angle=open_angle)

    return plot!(ps,xdata,ydata,color=:black,linewidth=1)
end

function add_arrow!(ps::Plots.Plot{Plots.PyPlotBackend},seg::Array{Float64,2},elements;num_arrows=1)

    if size(seg,1) <= 1
        return ps
    end

    lastrow = size(seg,1)
    if seg[end,:] == seg[end-1,:]
        lastrow -= 1
    end
    firstrow = 1
    if lastrow-firstrow >= 2
      firstrow += 1
    elseif lastrow-firstrow <= 0
      return ps
    end
    deg = (lastrow-firstrow) < 3 ? 1 : 3

    spl = ParametricSpline(seg[firstrow:lastrow,:]',k=deg)

    interv = 1/(num_arrows+1)
    for f in interv:interv:1-interv
        xhead, yhead = spl(f)
        w = induce_velocity(xhead+im*yhead,elements,0)
        arrowhead!(ps,xhead,yhead,real(w),imag(w),0.1)
    end
    return ps
end

"""
    add_arrow!(p::Plot,level::Integer,elements,[num_arrows=2])

Add arrow(s) on a given contour level of the plot `p`. The `elements` are the potential flow element(s)
that generate the flow. The number of arrows can be specified with the optional argument, but defaults
to 2
"""
function add_arrow!(ps::Plots.Plot{Plots.PyPlotBackend},level::Integer,el;a...)
    seg_list = get_segment_coords(ps,level)
    if length(seg_list) > 0
        for seg in seg_list
            add_arrow!(ps,seg,el;a...)
        end
    else
        return ps
    end
end

function add_arrows!(ps::Plots.Plot{Plots.PyPlotBackend},el;a...)
    contourset = get_series_handle(ps::Plots.Plot)
    nlev = length(contourset."levels")
    for lev in 1:2:nlev
        add_arrow!(ps,lev,el;a...)
    end
    return ps
end

####

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
      plan_constraints(u,t) = ViscousFlow.TimeMarching.plan_constraints(u,t,sys)
      r₁(u,t) = ViscousFlow.TimeMarching.r₁(u,t,sys)
      r₂(u,t) = ViscousFlow.TimeMarching.r₂(u,t,sys)

      return IFHERK(state,f,sys.Δt,plan_intfact,plan_constraints,(r₁,r₂),rk=ViscousFlow.TimeMarching.RK31,isstored=true), sys, state, f
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
     compute_trajectory(u,v,X₀::Vector,Tmax,Δt)

  Calculate the trajectory of a tracer particle with initial location(s) `X₀`, which
  can be specified as either a single vector `[x0,y0]` or a vector of vectors
  for multiple tracer particles. The arguments
  `u` and `v` are interpolated velocity field components, `Tmax` is the final
  integration time, and `Δt` is the time step size. The output is the solution
  structure for the `OrdinaryDiffEq` package (or, for multiple particles, a vector
  of such solution structures).
  """
  function compute_trajectory(ufield::AbstractInterpolation{T,2},
                              vfield::AbstractInterpolation{T,2},
                              X₀::Vector{S},Tmax::Real,Δt::Real) where {T,S<:Real}

    vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

    sol = _solve_trajectory(vfcn!,X₀,Tmax,Δt)
    return sol

  end

  function compute_trajectory(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},
     pts::Vector{Vector{S}},Tmax,Δt) where {T,S<:Real}

    sol_array = ODESolution[]
    for X₀ in pts
      sol = compute_trajectory(ufield,vfield,X₀,Tmax,Δt)
      push!(sol_array,sol)
    end
    return sol_array

  end


  """
     compute_trajectory(vel::Edges,sys,X₀::Vector/Vector{Vector},Tmax,Δt)

  Calculate the trajectory of a particle with initial location `X₀`. The argument
  `vel` is edge-type grid data, `sys` is a Navier-Stokes type system, `Tmax` is the final
  integration time, and `Δt` is the time step size. The output is the solution
  structure for the `OrdinaryDiffEq` package.
  """
  compute_trajectory(u::Edges, sys::NavierStokes, X₀::Union{Vector{S},Vector{Vector{S}}},Tmax,Δt) where S <: Real =
      compute_trajectory(interpolatable_field(u,sys.grid)...,X₀,Tmax,Δt)

  """
      compute_trajectory(elements,X₀::Vector,Tmax,Δt)

  Calculate the trajectory of a particle with initial location `X₀`. The
  argument `elements` is a potential flow `Element` type or group of `Element` types.
  `Tmax` is the final integration time, and `Δt` is the time step size. The output is the solution
  structure for the `OrdinaryDiffEq` package.
  """
  function compute_trajectory(elements, X₀::Vector{S}, Tmax, Δt) where S <: Real

    function vfcn!(dR,R,p,t)
      dR_complex = induce_velocity(R[1]+im*R[2], elements, t)
      dR[1] = real(dR_complex)
      dR[2] = imag(dR_complex)
      return dR
    end

    sol = _solve_trajectory(vfcn!,X₀,Tmax,Δt)
    return sol

  end

  function _solve_trajectory(vfcn!,u0,Tmax,Δt)
    Path = ODEProblem(vfcn!,u0,(0.0,Tmax))
    sol = solve(Path,ABM54(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)
  end


 """
     compute_trajectories(elements,tracer_start,Tmax,Δt)

 Calculate the trajectories of a set of tracer particles in a potential flow
 generated by the elements `elements`. The tracers' initial positions are
 provided as a vector of complex positions in `tracer_start`. The final time
 and time step are provided as `Tmax` and `Δt`. The output is a tuple of arrays
 of the x and y coordinates of the trajectories. Each column of these
 arrays corresponds to a single tracer history.
 """
 function compute_trajectories(elements, tracer_start::Vector{<:Number}, Tmax, Δt)

     tracer_x = []
     tracer_y = []
     for z0 in tracer_start
         sol = compute_trajectory(elements,[real(z0),imag(z0)],Tmax,Δt)
         xy = transpose(hcat(sol.u...))
         push!(tracer_x,xy[:,1])
         push!(tracer_y,xy[:,2])
     end
     x = hcat(tracer_x...)
     y = hcat(tracer_y...)

     return x, y

 end



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
 function mag(u::Edges{C}) where {C <: ViscousFlow.CartesianGrids.CellType}

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


### Boundary layer routines

function fsrhs!(du,u,p,t)
  _, m, _ = p

  du[1] = u[2]
  du[2] = u[3]
  du[3] = -0.5*(m+1)*u[1]*u[3] - m*(1-u[2]^2)
end

function fs_integrate(h0,p)
  Vw, m, ηmax = p

  f0 = -2Vw
  g0 = 0.0
  gL = 1.0

  u0 = [f0;g0;h0]
  ηspan = (0.0,ηmax)


  prob = ODEProblem(fsrhs!,u0,ηspan,p)
  sol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-15)

  resid = sol[2,end] - gL

  return resid, sol

end

"""
    falknerskan(β[,Vw=0][,ηmax=8][,h0init=1.3]) -> u, η, d99, dstar, theta, cf

Compute the Falkner-Skan boundary layer solution for parameter `β`,
where `β` can be a value larger than -0.1999 (the separation case).
It returns, in order
* the velocity profile (\$u/U\$)
* the corresponding vertical coordinates \$\\eta\$
* the proportionality factor on the 99 percent thickness
* the proportionality factor on the displacement thickness
* the proportionality factor on the momentum thickness
* the proportionality factor on the skin friction coefficient
"""
function falknerskan(β;Vw = 0.0,ηmax = 8.0,h0init=1.3)

  m = β/(2-β)
  p = [Vw,m,ηmax]

  h0 = find_zero(x -> fs_integrate(x,p)[1],h0init)
  resid,sol = fs_integrate(h0,p)

  return _fs_velocity(sol), _fs_eta(sol), blthickness_99(sol), blthickness_displacement(sol),
  blthickness_momentum(sol), blskinfriction(sol)

end

_fs_streamfunction(sol::OrdinaryDiffEq.ODESolution) = sol[1,:]
_fs_velocity(sol::OrdinaryDiffEq.ODESolution) = sol[2,:]
_fs_eta(sol::OrdinaryDiffEq.ODESolution) = sol.t

function blthickness_99(sol::OrdinaryDiffEq.ODESolution)

    u, η = _fs_velocity(sol), _fs_eta(sol)

    sign_diff = sign.(u .- 0.99)  # +1 where u > 0.99, -1 otherwise
    i0 = findfirst(diff(sign_diff) .== 2.0)
    return sol.t[i0]+ (0.99-u[i0])/(u[i0+1]-u[i0])*(η[i0+1]-η[i0])
end

function blthickness_displacement(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum((1 .- uhalf).*diff(η))
end

function blthickness_momentum(sol::OrdinaryDiffEq.ODESolution)
    u, η = _fs_velocity(sol), _fs_eta(sol)

    uhalf = 0.5*(u[2:end].+u[1:end-1])
    return sum(uhalf.*(1 .- uhalf).*diff(η))
end

blskinfriction(sol::OrdinaryDiffEq.ODESolution) = 2*sol[3,1]

end
