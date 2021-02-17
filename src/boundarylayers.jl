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
function falknerskan(β;Vw = 0.0,ηmax = 10.0,h0init=1.3)

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
