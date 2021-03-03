####### RAYLEIGH FLOW ########


HeatFlux(h01::StagnationEnthalpy,h02::StagnationEnthalpy) = HeatFlux(h02-h01)

function HeatFlux(T01::StagnationTemperature,T02::StagnationTemperature;gas::PerfectGas=DefaultPerfectGas)
  h01 = StagnationEnthalpy(T01;gas=gas)
  h02 = StagnationEnthalpy(T02;gas=gas)
  return HeatFlux(h01,h02)
end

function T0OverT0Star(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return TemperatureRatio(((γ+1)*M^2*(2+(γ-1)*M^2))/(1+γ*M^2)^2)
end

function TOverTStar(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return TemperatureRatio(((γ+1)^2*M^2)/(1+γ*M^2)^2)
end

function MachNumber(T0_over_T0star::TemperatureRatio,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
   value(T0_over_T0star) <= 1.0 || error("T0/T0* must be 1 or smaller")

    Msub = find_zero(x -> T0OverT0Star(MachNumber(x),RayleighFlow,gas=gas)-T0_over_T0star,(0.001,1),order=16)
    max_T0ratio = T0OverT0Star(MachNumber(1e10),RayleighFlow,gas=gas)

    if value(T0_over_T0star) < value(max_T0ratio)
        return MachNumber(Msub)
    else
        Msup = find_zero(x -> T0OverT0Star(MachNumber(x),RayleighFlow,gas=gas)-T0_over_T0star,(1,1e10),order=16)
        return MachNumber(Msub), MachNumber(Msup)
    end
end

function POverPStar(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return PressureRatio((γ+1)/(1+γ*M^2))
end

function VOverVStar(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return VelocityRatio(((γ+1)*M^2)/(1+γ*M^2))
end

function ρOverρStar(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  return DensityRatio(1/VOverVStar(M,RayleighFlow,gas=gas))
end

function P0OverP0Star(M::MachNumber,::Type{RayleighFlow};gas::PerfectGas=DefaultPerfectGas)
  γ = SpecificHeatRatio(gas)
  return PressureRatio((γ+1)/(1+γ*M^2)*((2+(γ-1)*M^2)/(γ+1))^(γ/(γ-1)))
end
