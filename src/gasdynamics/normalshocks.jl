###### NORMAL SHOCK RELATIONS #######

function MachNumber(M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    M2sq = (1 + 0.5*(γ-1)*M1^2)/(γ*M1^2-0.5*(γ-1))
    return MachNumber(sqrt(M2sq))
end

function TemperatureRatio(M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    return TemperatureRatio(1 + 2*(γ-1)/(γ+1)^2*((1+γ*M1^2)/M1^2)*(M1^2-1))
end

function PressureRatio(M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    return PressureRatio(1 + 2*γ/(γ+1)*(M1^2-1))
end

function DensityRatio(M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    return DensityRatio((γ+1)*M1^2/(2+(γ-1)*M1^2))
end

function StagnationPressureRatio(M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    M2 = MachNumber(M1,NormalShock,gas=gas)
    pratio = PressureRatio(M1,NormalShock,gas=gas)
    p01_unit = StagnationPressure(Pressure(1),M1,Isentropic,gas=gas)
    p02_unit = StagnationPressure(Pressure(1),M2,Isentropic,gas=gas)
    return StagnationPressureRatio(pratio*p02_unit/p01_unit)
end

function Entropy(s1::Entropy,M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    cv = SpecificHeatVolume(gas)
    ds = cv*(log(1+2*γ/(γ+1)*(M1^2-1)) - γ*log((γ+1)*M1^2/(2+(γ-1)*M1^2)))
    return Entropy(s1+ds)

end
