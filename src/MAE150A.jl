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
  using Conda


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



      if isdefined(Main, :IJulia) && Main.IJulia.inited
        # The Pkg.build does not work if non-development package, so need to
        # ensure JIT install of matplotlib using Conda
        _hasmatplotlib() || Conda.add("matplotlib")
      else
        # For develop (e.g., CI), force re-build of PyCall with internal Python dist,
        # to make sure matplotlib is installed:
        if _isuserwritable(proj_file)
          ENV["PYTHON"] = ""
          Pkg.build("PyCall")
        else
          error("Project file is not writable. Cannot build PyCall")
        end
      end


      # Get LaTeXStrings from PyPlot
      #using PyPlot: LaTeXStrings
      using LaTeXStrings

      Plots.pyplot()
      rcParams = Plots.PyPlot.PyDict(Plots.PyPlot.matplotlib."rcParams")

      # Ensure that LaTeX stuff is handled
      rcParams["mathtext.fontset"] = "cm"

      Plots.default(markerstrokealpha = 0, legend = false,
        dpi = 100, size = (400, 300), grid = false)

      include("arrows.jl")

    end

  end

  _isuserwritable(file) = (uperm(file) >> 1) & 1 != 0
  _hasmatplotlib() = haskey(Conda._installed_packages_dict(),"matplotlib")

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


include("fileio.jl")
include("trajectories.jl")
include("boundarylayers.jl")
include("potentialflow.jl")


end
