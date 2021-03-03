######## ISENTROPIC RELATIONS #########

function TemperatureRatio(pratio::PressureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    TemperatureRatio(pratio^((γ-1)/γ))
end

function TemperatureRatio(ρratio::DensityRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    TemperatureRatio(ρratio^(γ-1))
end

function PressureRatio(Tratio::TemperatureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    PressureRatio(Tratio^(γ/(γ-1)))
end

function P0OverP(M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    Tratio = T0OverT(M,Isentropic,gas=gas)
    PressureRatio(Tratio,Isentropic,gas=gas)
end


function DensityRatio(Tratio::TemperatureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    DensityRatio(Tratio^(1/(γ-1)))
end

function ρ0Overρ(M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    Tratio = T0OverT(M,Isentropic,gas=gas)
    DensityRatio(Tratio,Isentropic,gas=gas)
end

# Temperature relations
function StagnationTemperature(T::Temperature,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    StagnationTemperature(T*T0OverT(M,Isentropic,gas=gas))
end

StagnationTemperature(T::Temperature,M::MachNumber;gas::PerfectGas=DefaultPerfectGas) = StagnationTemperature(T,M,Isentropic,gas=gas)

Temperature(T0::StagnationTemperature,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas) = Temperature(T0/T0OverT(M,gas=gas))

Temperature(T0::StagnationTemperature,M::MachNumber;gas::PerfectGas=DefaultPerfectGas) = Temperature(T0,M,Isentropic,gas=gas)

function MachNumber(T_over_T0::TemperatureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    value(T_over_T0) <= 1.0 || error("T/T0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    M2 = ((1.0/T_over_T0)-1)*2/(γ-1)
    MachNumber(sqrt(M2))
end

MachNumber(T_over_T0::TemperatureRatio;gas::PerfectGas=DefaultPerfectGas) = MachNumber(T_over_T0,Isentropic,gas=gas)

MachNumber(T::Temperature,T0::StagnationTemperature,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas) = MachNumber(TemperatureRatio(T/T0),gas=gas)

MachNumber(T::Temperature,T0::StagnationTemperature;gas::PerfectGas=DefaultPerfectGas) = MachNumber(T,T0,Isentropic,gas=gas)

# Pressure relations
function StagnationPressure(p::Pressure,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    StagnationPressure(p*PressureRatio(T0OverT(M,gas=gas),Isentropic,gas=gas))
end
function Pressure(p0::StagnationPressure,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    Pressure(p0/PressureRatio(T0OverT(M,gas=gas),Isentropic,gas=gas))
end

function MachNumber(p_over_p0::PressureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    value(p_over_p0) <= 1.0 || error("p/p0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    MachNumber(TemperatureRatio(p_over_p0,Isentropic,gas=gas))
end

MachNumber(p::Pressure,p0::StagnationPressure,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas) = MachNumber(PressureRatio(p/p0),gas=gas)

# Density relations
function StagnationDensity(ρ::Density,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    StagnationDensity(ρ*DensityRatio(T0OverT(M,gas=gas),Isentropic,gas=gas))
end
function Density(ρ0::StagnationDensity,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    Density(ρ0/DensityRatio(T0OverT(M,gas=gas),Isentropic,gas=gas))
end

function MachNumber(ρ_over_ρ0::DensityRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    value(ρ_over_ρ0) <= 1.0 || error("ρ/ρ0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    MachNumber(TemperatureRatio(ρ_over_ρ0,Isentropic,gas=gas))
end

MachNumber(ρ::Density,ρ0::StagnationDensity,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas) = MachNumber(DensityRatio(ρ/ρ0),gas=gas)


# Area-Mach number relation
function AOverAStar(M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    return AreaRatio(1/M*(2/(γ+1)*(T0OverT(M,gas=gas)))^(0.5*(γ+1)/(γ-1)))
end
AOverAStar(M::MachNumber;gas::PerfectGas=DefaultPerfectGas) = AOverAStar(M,Isentropic,gas=gas)

function AStar(A::Area,M::MachNumber,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    return Area(A/AOverAStar(M,gas=gas))
end
AStar(A::Area,M::MachNumber;gas::PerfectGas=DefaultPerfectGas) = AStar(A,M,Isentropic,gas=gas)

function SubsonicMachNumber(A_over_Astar::AreaRatio,::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas)
    value(A_over_Astar) >= 1.0 || error("A/A* must be 1 or larger")
    Msub = find_zero(x -> AOverAStar(MachNumber(x))-A_over_Astar,(0,1),order=16)
    return MachNumber(Msub)
end
function SupersonicMachNumber(A_over_Astar::AreaRatio,::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas)
  value(A_over_Astar) >= 1.0 || error("A/A* must be 1 or larger")
  Msup = find_zero(x -> AOverAStar(MachNumber(x))-A_over_Astar,(1,Inf),order=16)
  return MachNumber(Msup)
end

function SubsonicMachNumber(A::Area,Aloc::Area,ploc_over_p0::PressureRatio,::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas)
    Mloc = MachNumber(ploc_over_p0,Isentropic,gas=gas)
    Astar = AStar(Aloc,Mloc,gas=gas)
    SubsonicMachNumber(A,Astar,Isentropic,gas=gas)
end

function SupersonicMachNumber(A::Area,Aloc::Area,ploc_over_p0::PressureRatio,::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas)
    Mloc = MachNumber(ploc_over_p0,Isentropic,gas=gas)
    Astar = AStar(Aloc,Mloc,gas=gas)
    SupersonicMachNumber(A,Astar,Isentropic,gas=gas)
end

SubsonicMachNumber(A::Area,Astar::Area, ::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas) =
    SubsonicMachNumber(AreaRatio(A/Astar),Isentropic,gas=gas)

SupersonicMachNumber(A::Area,Astar::Area, ::Type{Isentropic}; gas::PerfectGas=DefaultPerfectGas) =
    SupersonicMachNumber(AreaRatio(A/Astar),Isentropic,gas=gas)


function MachNumber(A_over_Astar::AreaRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    return SubsonicMachNumber(A_over_Astar,Isentropic,gas=gas),
           SupersonicMachNumber(A_over_Astar,Isentropic,gas=gas)
end

MachNumber(A_over_Astar::AreaRatio;gas::PerfectGas=DefaultPerfectGas) = MachNumber(A_over_Astar,Isentropic,gas=gas)

function MachNumber(M1::MachNumber,A1::Area,A2::Area,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    M2sub, M2sup = MachNumber(AreaRatio(A2/AStar(A1,M1,gas=gas)),gas=gas)
end
MachNumber(M1::MachNumber,A1::Area,A2::Area;gas::PerfectGas=DefaultPerfectGas) = MachNumber(M1,A2,A2,Isentropic,gas=gas)


function SubsonicPOverP0(A::Area,Aloc::Area,ploc_over_p0::PressureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
  Mloc = MachNumber(ploc_over_p0,Isentropic,gas=gas)
  Astar = AStar(Aloc,Mloc,gas=gas)
  SubsonicPOverP0(A,Astar,Isentropic,gas=gas)
end

function SupersonicPOverP0(A::Area,Aloc::Area,ploc_over_p0::PressureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
  Mloc = MachNumber(ploc_over_p0,Isentropic,gas=gas)
  Astar = AStar(Aloc,Mloc,gas=gas)
  SupersonicPOverP0(A,Astar,Isentropic,gas=gas)
end

function SubsonicPOverP0(A::Area,Astar::Area,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
  M = SubsonicMachNumber(A,Astar,Isentropic,gas=gas)
  PressureRatio(1/P0OverP(M,Isentropic,gas=gas))
end

function SupersonicPOverP0(A::Area,Astar::Area,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
  M = SupersonicMachNumber(A,Astar,Isentropic,gas=gas)
  PressureRatio(1/P0OverP(M,Isentropic,gas=gas))
end
