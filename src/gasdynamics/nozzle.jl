# Converging-diverging nozzles
export Nozzle, areas, positions, converging, diverging, throat,
      pressure_nozzle, machnumber_nozzle, nozzle_quality


struct Nozzle{T,S}
    Ai :: T
    At :: T
    Ae :: T
    x :: S
    areas :: Vector{T}
    throat_index :: Integer
end

function Nozzle(Ai::Area,At::Area,Ae::Area;xmin=-1.0,xmax=3.0,len=201)
    x = range(xmin,xmax,length=len)
    a = areas(x,Ai,At,Ae)
    Nozzle(Ai,At,Ae,x,a,_find_throat(a,At))
end

function areas(x,Ai::Area,At::Area,Ae::Area)
    a = Array{Area,1}(undef,length(x))
    araw = _nozzle_area(x,value(Ai),value(At),value(Ae))
    a .= Area.(araw)
end

function _find_throat(a::AbstractVector{T},At::Area) where T <: Area
    indx = findall(x -> (value(x) ≈ value(At)),a)
    isempty(indx) && error("Can't find the throat")
    return indx[1]
end

areas(n::Nozzle) = n.areas
positions(n::Nozzle) = n.x
throat(n::Nozzle) = n.At
exit(n::Nozzle) = n.Ae
inlet(n::Nozzle) = n.Ai
converging(n::Nozzle) = view(areas(n),1:n.throat_index)
diverging(n::Nozzle) = view(areas(n),n.throat_index+1:length(areas(n)))

function Base.show(io::IO, noz::Nozzle)
  println(io, "Converging-diverging nozzle")
  println(io, "   Inlet area (sq cm)= "*string(value(inlet(noz),SqCM)))
  println(io, "   Throat area (sq cm)= "*string(value(throat(noz),SqCM)))
  println(io, "   Exit area (sq cm)= "*string(value(exit(noz),SqCM)))
end

# Shaping of nozzle cross-section
_rad_conv(x,Rt,Ri) = Rt + (Ri-Rt)*x^2
_rad_div(x,Rt,Re,c) = Re + (Rt-Re)*exp(-c*x^2)
_radius_shape(x,Ri,Rt,Re,c) = x <= 0 ? _rad_conv(x,Rt,Ri) : _rad_div(x,Rt,Re,c)
function _nozzle_area(x,Ai,At,Ae;c=1)
    return π*_radius_shape.(x,sqrt(Ai/π),sqrt(At/π),sqrt(Ae/π),c).^2
end

# Diverging nozzle processes

function pressure_diverging_nozzle(A::AbstractVector{T},Astar::Area,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    p = Pressure[]
    p_over_p0 = SupersonicPOverP0(A,Astar,Isentropic,gas=gas)
    for pr in p_over_p0
        push!(p,Pressure(pr*p0))
    end

    return p
end

function machnumber_diverging_nozzle(A::AbstractVector{T},Astar::Area;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    return SupersonicMachNumber(A,Astar,Isentropic,gas=gas)
end


function pressure_diverging_nozzle_with_shock(A::AbstractVector{T},As::Area,Astar1::Area,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    i1 = _index_upstream_of_shock(A,As)

    p = Pressure[]
    pup_over_p01 = SupersonicPOverP0(A[1:i1],Astar1,Isentropic,gas=gas)
    for pr in pup_over_p01
        push!(p,Pressure(pr*p01))
    end

    M1, M2, Astar2, p02_over_p01 = _shock(As,Astar1)

    p02 = StagnationPressure(p02_over_p01*p01)

    pdown_over_p02 = SubsonicPOverP0(A[i1+1:end],Astar2,Isentropic,gas=gas)
    for pr in pdown_over_p02
        push!(p,Pressure(pr*p02))
    end
    return p
end

function machnumber_diverging_nozzle_with_shock(A::AbstractVector{T},As::Area,Astar1::Area,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    i1 = _index_upstream_of_shock(A,As)

    M = SupersonicMachNumber(A[1:i1],Astar1,Isentropic,gas=gas)
    M1, M2, Astar2, _ = _shock(As,Astar1)

    Mdown = SubsonicMachNumber(A[i1+1:end],Astar2,Isentropic,gas=gas)
    for Mi in Mdown
        push!(M,Mi)
    end
    return M
end

# Given an exit area Ae, shock area As, and Astar1, p01 in upstream nozzle
# return the pressure pe at exit
function pressure_diverging_nozzle_with_shock(Ae::Area,As::Area,Astar1::Area,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    value(Astar1) <= value(As) <= value(Ae) || error("Shock area outside of bounds")

    M1, M2, Astar2, p02_over_p01 = _shock(As,Astar1)

    # conditions downstream of shock
    p02 = StagnationPressure(p02_over_p01*p01)

    pe_over_p02 = SubsonicPOverP0(Ae,Astar2,Isentropic,gas=gas)
    return Pressure(pe_over_p02*p02)
end

# Given an exit area Ae, shock area As, and Astar1, p01 in upstream nozzle
# return the Mach number Me at exit
function machnumber_diverging_nozzle_with_shock(Ae::Area,As::Area,Astar1::Area,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    value(Astar1) <= value(As) <= value(Ae) || error("Shock area outside of bounds")

    M1, M2, Astar2, _ = _shock(As,Astar1)

    return SubsonicMachNumber(Ae,Astar2,Isentropic,gas=gas)
end

function _shock(As::Area,Astar1::Area;gas::PerfectGas=DefaultPerfectGas)
    M1 = SupersonicMachNumber(AreaRatio(As/Astar1),Isentropic,gas=gas)
    M2 = MachNumber(M1,NormalShock,gas=gas)
    Astar2 = AStar(As,M2,Isentropic,gas=gas)
    p02_over_p01 = StagnationPressureRatio(M1,NormalShock,gas=gas)
    return M1, M2, Astar2, p02_over_p01
end


# Find the index in array A corresponding to just before As.
# Assumes that A is diverging nozzle.
function _index_upstream_of_shock(A::AbstractVector{T},As::Area) where T <: Area
    value(A[1]) <= value(As) <= value(A[end]) || error("Shock area outside of bounds")
    indx = findall(x -> x >= 1.0,diff(sign.(value.(A) .- value(As))))
    return indx[end]
end

function machnumber_nozzle_with_shock(noz::Nozzle,pb::Pressure,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    As = _shock_area(noz,pb,p01,gas=gas)
    vcat(machnumber_subsonic(converging(noz),converging(noz)[end],gas=gas),
         machnumber_diverging_nozzle_with_shock(diverging(noz),As,throat(noz),p01,gas=gas))
end

function pressure_nozzle_with_shock(noz::Nozzle,pb::Pressure,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    As = _shock_area(noz,pb,p01,gas=gas)
    vcat(pressure_subsonic(converging(noz),converging(noz)[end],p01,gas=gas),
         pressure_diverging_nozzle_with_shock(diverging(noz),As,throat(noz),p01,gas=gas))
end

function machnumber_nozzle_supersonic(noz::Nozzle;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    vcat(machnumber_subsonic(converging(noz),converging(noz)[end],gas=gas),
         machnumber_diverging_nozzle(diverging(noz),throat(noz),gas=gas))
end

function pressure_nozzle_supersonic(noz::Nozzle,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
    vcat(pressure_subsonic(converging(noz),converging(noz)[end],p0,gas=gas),
         pressure_diverging_nozzle(diverging(noz),throat(noz),p0,gas=gas))
end

function machnumber_subsonic(A::AbstractVector{T},pb::Pressure,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
  Mb = MachNumber(PressureRatio(pb/p0),Isentropic,gas=gas)
  Astar = AStar(A[end],Mb,gas=gas)
  machnumber_subsonic(A,Astar,gas=gas)
end

function machnumber_subsonic(A::AbstractVector{T},Astar::Area;gas::PerfectGas=DefaultPerfectGas) where T <: Area
  return SubsonicMachNumber(A,Astar,Isentropic,gas=gas)
end

function pressure_subsonic(A::AbstractVector{T},pb::Pressure,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
  Mb = MachNumber(PressureRatio(pb/p0),Isentropic,gas=gas)
  Astar = AStar(A[end],Mb,gas=gas)
  pressure_subsonic(A,Astar,p0,gas=gas)
end

function pressure_subsonic(A::AbstractVector{T},Astar::Area,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) where T <: Area
  p_over_p0 = SubsonicPOverP0(A,Astar,Isentropic,gas=gas)
  p = Pressure[]
  for pr in p_over_p0
      push!(p,Pressure(pr*p0))
  end
  return p
end

function machnumber_nozzle(noz::Nozzle,pb::Pressure,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas)
  pb_subcrit = _pressure_subsonic_critical(noz,p01,gas=gas)
  pb_supcrit = _pressure_supersonic_shock_critical(noz,p01,gas=gas)

  M = pb > pb_subcrit ? machnumber_subsonic(areas(noz),pb,p01,gas=gas) :
            (pb > pb_supcrit ? machnumber_nozzle_with_shock(noz,pb,p01,gas=gas) : machnumber_nozzle_supersonic(noz,gas=gas))

  return M
end

function pressure_nozzle(noz::Nozzle,pb::Pressure,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas)
  pb_subcrit = _pressure_subsonic_critical(noz,p01,gas=gas)
  pb_supcrit = _pressure_supersonic_shock_critical(noz,p01,gas=gas)
  p = pb > pb_subcrit ? pressure_subsonic(areas(noz),pb,p01,gas=gas) :
            (pb > pb_supcrit ? pressure_nozzle_with_shock(noz,pb,p01,gas=gas) : pressure_nozzle_supersonic(noz,p01,gas=gas))

  # Set the final pressure to the back pressure
  p[end] = pb

  return p
end



# Find area at which shock occurs
function _shock_area(noz::Nozzle,pb::Pressure,p01::StagnationPressure;gas::PerfectGas=DefaultPerfectGas)
  pb_crit = _pressure_subsonic_critical(noz,p01,gas=gas)
  pb < pb_crit || error("Back pressure too large for a shock")

  return Area(find_zero(x -> value(pressure_diverging_nozzle_with_shock(exit(noz),Area(x),throat(noz),p01))-value(pb),
    (value(diverging(noz)[1]),value(diverging(noz)[end-1])),order=8))
end

function nozzle_quality(noz::Nozzle,pb::Pressure,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas)
  pb_subcrit = _pressure_subsonic_critical(noz,p0,gas=gas)
  pb_supcrit_exitshock = _pressure_supersonic_shock_critical(noz,p0,gas=gas)
  pb_supcrit = _pressure_supersonic_critical(noz,p0,gas=gas)
  if pb > p0 || pb < 0.0
    return "Invalid back pressure"
  elseif pb > pb_subcrit
    return "Unchoked subsonic"
  elseif pb ≈ pb_subcrit
    return "Choked subsonic"
  elseif pb > pb_supcrit_exitshock
    return "Supersonic with normal shock"
  elseif pb > pb_supcrit
    return "Overexpanded supersonic"
  elseif pb ≈ pb_supcrit
    return "Perfectly expanded supersonic"
  else
    return "Underexpanded supersonic"
  end
  #println(pb_supcrit)
end

#=
function nozzle_quality(noz::Nozzle,pb::Pressure,p0::StagnationPressure;gas::PerfectGas=DefaultGasConstant)
    pb_subcrit = _pressure_subsonic_critical(noz,p0,gas=gas)
    pb_supcrit_exitshock = _pressure_supersonic_shock_critical(noz,p0,gas=gas)
    pb_supcrit = _pressure_supersonic_critical(noz,p0,gas=gas)
    nothing

    if pb > p0 || pb < 0.0
      return "Invalid back pressure"
    elseif pb > pb_subcrit
      return "Unchoked subsonic"
    elseif pb ≈ pb_subcrit
      return "Choked subsonic"
    elseif pb > pb_supcrit_exitshock
      return "Choked supersonic with normal shock"
    elseif pb > pb_supcrit
      return "Overexpanded supersonic"
    elseif pb ≈ pb_supcrit
      return "Perfectly expanded supersonic"
    else
      return "Underexpanded supersonic"
    end

end
=#

# Find critical back pressure for subsonic isentropic branch
function _pressure_subsonic_critical(noz::Nozzle,p0::StagnationPressure;gas::PerfectGas=DefaultGasConstant)
  pb_crit_over_p0 = SubsonicPOverP0(exit(noz),throat(noz),Isentropic,gas=gas)
  return Pressure(pb_crit_over_p0*p0)
end

# Find critical back pressure for supersonic isentropic branch with shock at exit
function _pressure_supersonic_shock_critical(noz::Nozzle,p0::StagnationPressure;gas::PerfectGas=DefaultGasConstant)
  pe_supcrit = _pressure_supersonic_critical(noz,p0,gas=gas)
  Me_supcrit = SupersonicMachNumber(exit(noz),throat(noz),Isentropic,gas=gas)

  pb_over_pe_supcrit = PressureRatio(Me_supcrit,NormalShock,gas=gas)
  return Pressure(pb_over_pe_supcrit*pe_supcrit)
end

# Find critical back pressure for supersonic isentropic branch
function _pressure_supersonic_critical(noz::Nozzle,p0::StagnationPressure;gas::PerfectGas=DefaultGasConstant)
  pe_supcrit_over_p0 = SupersonicPOverP0(exit(noz),throat(noz),Isentropic,gas=gas)
  return Pressure(pe_supcrit_over_p0*p0)
end
