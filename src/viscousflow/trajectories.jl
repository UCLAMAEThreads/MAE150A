
"""
   compute_trajectory(vel::Edges,sys,X₀::Vector/Vector{Vector},Tmax[,Δt=0.001])

Calculate the trajectory of a particle with initial location `X₀`. The argument
`vel` is edge-type grid data, `sys` is a Navier-Stokes type system, `Tmax` is the final
integration time, and `Δt` is the time step size. The output is the solution
structure for the `OrdinaryDiffEq` package.
"""
compute_trajectory(u::ViscousFlow.Edges, sys::ViscousFlow.ILMSystem, X₀::Union{Vector{S},Vector{Vector{S}}},Tmax;Δt=DEFAULT_DT) where S <: Real =
    compute_trajectory(ViscousFlow.interpolatable_field(u,sys.base_cache.g)...,X₀,Tmax;Δt=Δt)


"""
  field_along_trajectory(f::GridData,sys::NavierStokes,traj::ODESolution[,deriv=0])

Evaluate field `f` (given as grid data) along the trajectory specified by `traj`.
The output is the history of `f` along this trajectory. If `f` is a vector field,
then the component histories are output as a tuple. If `deriv=1`, then it
computes the time derivative of the field along the trajectory. The default
is `deriv=0` (no derivative).
"""
field_along_trajectory(d::ViscousFlow.GridData,sys,traj;deriv=0) = _field_along_trajectory(d,sys,traj,Val(deriv))

function _field_along_trajectory(v::ViscousFlow.VectorGridData,sys::ViscousFlow.ILMSystem,traj::ODESolution,::Val{0})
 vfield_x, vfield_y = ViscousFlow.interpolatable_field(v,sys.base_cache.g)

 vx_traj = eltype(v)[]
 vy_traj = eltype(v)[]
 for x in traj.u
   push!(vx_traj,vfield_x(x...))
   push!(vy_traj,vfield_y(x...))
 end

 return vx_traj, vy_traj
end

function _field_along_trajectory(s::ViscousFlow.ScalarGridData,sys::ViscousFlow.ILMSystem,traj::ODESolution,::Val{0})
 sfield = ViscousFlow.interpolatable_field(s,sys.base_cache.g)

 s_traj = eltype(sfield)[]
 for x in traj.u
   push!(s_traj,sfield(x...))
 end

 return s_traj
end

function _field_along_trajectory(v::ViscousFlow.VectorGridData,sys::ViscousFlow.ILMSystem,traj::ODESolution,::Val{1})
    utraj, vtraj = _field_along_trajectory(v,sys,traj,Val(0))
    return ddt(utraj,traj.t), ddt(vtraj,traj.t)
end

_field_along_trajectory(s::ViscousFlow.ScalarGridData,sys::ViscousFlow.ILMSystem,traj::ODESolution,::Val{1}) =
    ddt(_field_along_trajectory(s,sys,traj,Val(0)),traj.t)
