module MAE150A

  using Reexport

  @reexport using ViscousFlow
  @reexport using Plots
  @reexport using OrdinaryDiffEq


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

end
