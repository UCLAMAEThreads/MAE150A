"""
    save_ns_solution(filen::String,integrator)

Save a state of the solution in a file with the provided
filename `filen`.
"""
function save_ns_solution(filen,integrator)

    sys = integrator.p
    save(filen,"u",integrator.u,
               "t",integrator.t,
               "prob",ViscousFlow.ImmersedLayers.regenerate_problem(sys,sys.base_cache.bl))
    return nothing

end

"""
      load_ns_solution(filen::String)

  Given a JLD file with name `filen`, load in the Navier-Stokes solution data
  stored in this file and set up various solution variables. An example:

  ```
  u, t, sys = load_ns_solution("myfile.jld2")
  ```

  In this example, `u` is the flow state vector, `t` the time, `sys` is the NS system
  metadata and operators.
"""
function load_ns_solution(filen)

  d = load(filen)
  u = d["u"]
  t = d["t"]
  sys = ViscousFlow.ImmersedLayers.construct_system(d["prob"])

  return u, t, sys
end
