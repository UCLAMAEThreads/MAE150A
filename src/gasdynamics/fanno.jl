######## FANNO FLOW ########


function FLStarOverD(M::MachNumber,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return FLOverD((1-M^2)/(γ*M^2) + (γ+1)/(2γ)*log((γ+1)*M^2/(2+(γ-1)*M^2)))
end

function MachNumber(fL_over_D::FLOverD,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
    Msub = find_zero(x -> FLStarOverD(MachNumber(x),FannoFlow,gas=gas)-fL_over_D,(0.001,1),order=16)
    max_fLD = FLStarOverD(MachNumber(1e10),FannoFlow,gas=gas)

    if value(fL_over_D) > value(max_fLD)
        return MachNumber(Msub)
    else
        Msup = find_zero(x -> FLStarOverD(MachNumber(x),FannoFlow,gas=gas)-fL_over_D,(1,1e10),order=16)
        return MachNumber(Msub), MachNumber(Msup)
    end
end

FLOverD(f::FrictionFactor,L::Length,D::Diameter) = FLOverD(f*L/D)
Length(fLoverD::FLOverD,D::Diameter,f::FrictionFactor) = Length(fLoverD*D/f)





function POverPStar(M::MachNumber,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return PressureRatio(1/M*sqrt((1+γ)/(2+(γ-1)*M^2)))
end

function ρOverρStar(M::MachNumber,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return DensityRatio(1/M*sqrt((2+(γ-1)*M^2)/(γ+1)))
end

function TOverTStar(M::MachNumber,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return TemperatureRatio((γ+1)/(2+(γ-1)*M^2))
end

function MachNumber(T_over_Tstar::TemperatureRatio,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
    M = find_zero(x -> TOverTStar(MachNumber(x),FannoFlow,gas=gas)-T_over_Tstar,(1e-10,1e10),order=16)
    return MachNumber(M)

end


function P0OverP0Star(M::MachNumber,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return StagnationPressureRatio(1/M*((2+(γ-1)*M^2)/(γ+1))^((γ+1)/(2(γ-1))))
end

function MachNumber(p0_over_p0star::StagnationPressureRatio,::Type{FannoFlow};gas::PerfectGas=DefaultPerfectGas)
    Msub = find_zero(x -> P0OverP0Star(MachNumber(x),FannoFlow,gas=gas)-p0_over_p0star,(0.001,1),order=16)
    max_p0_over_p0star = P0OverP0Star(MachNumber(1e10),FannoFlow,gas=gas)

    if value(p0_over_p0star) > value(max_p0_over_p0star)
        return MachNumber(Msub)
    else
        Msup = find_zero(x -> P0OverP0Star(MachNumber(x),FannoFlow,gas=gas)-p0_over_p0star,(1,1e10),order=16)
        return MachNumber(Msub), MachNumber(Msup)
    end
end
