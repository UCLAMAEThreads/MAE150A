# TRAJECTORY CALCULATION #

const DEFAULT_DT = 0.001

"""
   compute_trajectory(u,v,X₀::Vector,Tmax[,Δt=0.001])

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
                            X₀::Vector{S},Tmax::Real;Δt::Real=DEFAULT_DT) where {T,S<:Real}

  vfcn!(dR,R,p,t) = _vfcn_autonomous!(dR,R,p,t,ufield,vfield)

  sol = _solve_trajectory(vfcn!,X₀,Tmax,Δt)
  return sol

end

function compute_trajectory(ufield::AbstractInterpolation{T,2},vfield::AbstractInterpolation{T,2},
   pts::Vector{Vector{S}},Tmax;Δt=DEFAULT_DT) where {T,S<:Real}

  sol_array = ODESolution[]
  for X₀ in pts
    sol = compute_trajectory(ufield,vfield,X₀,Tmax,Δt=Δt)
    push!(sol_array,sol)
  end
  return sol_array

end


function _solve_trajectory(vfcn!,u0,Tmax,Δt)
  Path = ODEProblem(vfcn!,u0,(0.0,Tmax))
  sol = solve(Path,Tsit5(), dt = Δt, maxiters = 1e8, adaptive = false, dense = false)
end


function _vfcn_autonomous!(dR,R,p,t,u,v)
 dR[1] = u(R[1],R[2])
 dR[2] = v(R[1],R[2])

 return dR
end
