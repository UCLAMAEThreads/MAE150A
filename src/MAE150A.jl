module MAE150A

  using Pkg, InteractiveUtils, IJulia

  using Reexport

  @reexport using ViscousFlow
  @reexport using PotentialFlow
  @reexport using OrdinaryDiffEq
  #@reexport using LaTeXStrings

  #import Plots: plot

  using Interpolations
  using JLD
  using Requires
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

  function __init__()

    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin

      #ENV["PYTHON"] = ""
      #Pkg.add("PyCall")
      #Pkg.build("PyCall")
      #
      #Pkg.add("PyPlot")
      using PyPlot: PyCall, LaTeXStrings


      Plots.pyplot()
      rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

      # Ensure that LaTeX stuff is handled
      rcParams["mathtext.fontset"] = "cm"

      Plots.default(markerstrokealpha = 0, legend = false,
        dpi = 100, size = (400, 300), grid = false)

      include("arrows.jl")

    end

  end


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
    Base.eval(Main, Meta.parse("import IJulia"))
    path = joinpath(repo_directory,"notebook")
    IJulia.notebook(;dir=path)
  end




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
