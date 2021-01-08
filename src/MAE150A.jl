module MAE150A

  using Pkg, InteractiveUtils, IJulia

  using Reexport

  @reexport using ViscousFlow
  @reexport using PotentialFlow
  #@reexport using Plots
  @reexport using OrdinaryDiffEq
  @reexport using LaTeXStrings

  #import Plots: plot

  using Interpolations
  using JLD
  @reexport using RecursiveArrayTools
  using Dierckx
  using Roots
  #using PyCall
  #using PyPlot

  export initialize_environment,initialize_ns_solver,
        save_ns_solution,load_ns_solution, get_flowfield,
        compute_trajectory, compute_trajectories,
        field_along_trajectory, field_deriv_along_trajectory,
        convective_acceleration, mag, ddt, pressure,
        OseenVortex,
        complexgrid, vortex_patch,
        add_arrow!,add_arrows!,
        falknerskan


  repo_directory = joinpath(@__DIR__,"..")

  include("plot_recipes.jl")

  function tutorial_footer(; remove_homedir=true)
      display("text/markdown", """
      ## Appendix
       This lesson is part of the MAE150.jl repository, found at: <https://github.com/jdeldre/MAE150A>.
      """)
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
    #initialize_environment()
    Base.eval(Main, Meta.parse("import IJulia"))
    path = joinpath(repo_directory,"notebook")
    IJulia.notebook(;dir=path)
  end


  function initialize_environment()

    # Set the back end for Plots
    #pyplot()

    #rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

    # Ensure that LaTeX stuff is handled
    #rcParams["mathtext.fontset"] = "cm"

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
#=
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
=#

####


  """
      save_ns_solution(filen::String,integrator)

  Save a state of the Navier-Stokes solution in a file with the provided
  filename `filen`.
  """
  function save_ns_solution(filen,integrator)

      sys = integrator.p
      save(filen,"Re",sys.Re,"freestream",sys.U∞,"Δt",sys.Δt,"grid",sys.grid,
                "bodies",sys.bodies,"motions",sys.motions,
                 "u",integrator.u,"t",integrator.t,"motiontype",motiontype(sys))
      return nothing

  end

  motiontype(sys::NavierStokes{NX, NY, N, MT}) where {NX,NY,N,MT} = MT

  """
        load_ns_solution(filen::String)

    Given a JLD file with name `filen`, load in the Navier-Stokes solution data
    stored in this file and set up various solution variables. An example:

    ```
    u, t, sys = load_ns_solution("myfile.jld")
    ```

    In this example, `u` is the flow state vector, `t` the time, `sys` is the NS system
    metadata and operators.
  """
  function load_ns_solution(filen)

    d = load(filen)

    bodies = d["bodies"]
    motions = d["motions"]
    g = d["grid"]
    Δt = d["Δt"]
    U∞ = d["freestream"]
    Re = d["Re"]
    u = d["u"]
    t = d["t"]
    motiontype = d["motiontype"]

    sp = motiontype == ViscousFlow.StaticPoints ? true : false

    xlim = round.(limits(g,1),digits=15)
    ylim = round.(limits(g,2),digits=15)

    sys = NavierStokes(Re,cellsize(g),xlim,ylim,Δt,bodies,motions,
    freestream = U∞,static_points=sp)

    return u, t, sys
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
 function ddt(u::AbstractVector{T},t::AbstractVector{S};mydiff::Symbol=:forward_diff) where {T,S}
    du = eval(mydiff)(u)./eval(mydiff)(t)
    return du
 end

 # Some basic differencing routines
 function backward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end] .= u[2:end] .- u[1:end-1]
     du[1] = du[2]
     return du
 end
 function forward_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[1:end-1] .= u[2:end] .- u[1:end-1]
     du[end] = du[end-1]
     return du
 end
 function central_diff(u::AbstractVector{T}) where {T}
     du = zero(u)
     du[2:end-1] .= 0.5*u[3:end] .- 0.5*u[1:end-2]
     du[1] = du[2]
     du[end] = du[end-1]
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


include("trajectories.jl")
include("boundarylayers.jl")



end
