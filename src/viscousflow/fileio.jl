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

  sys = NavierStokes(Re,ViscousFlow.cellsize(g),xlim,ylim,Δt,bodies,motions,
  freestream = U∞,static_points=sp)

  return u, t, sys
end
