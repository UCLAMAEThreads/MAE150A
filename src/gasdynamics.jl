#=
Quasi-1d gas dynamics routines
=#

import Base:+,*,-,/,^

using Roots

######## UNITS #########
abstract type ThermodynamicUnits end
abstract type SI <: ThermodynamicUnits end
abstract type SIOtherUnits <: ThermodynamicUnits end
abstract type Imperial <: ThermodynamicUnits end
abstract type OtherUnits <: ThermodynamicUnits end

abstract type Dimensionless <: ThermodynamicUnits end

function default_unit end

convert_unit(::Type{U},::Type{U},val) where {U<:ThermodynamicUnits} = val

# temperature
abstract type Kelvin <: SI end
abstract type K <: SIOtherUnits end
abstract type Celsius <: SIOtherUnits end
abstract type C <: SIOtherUnits end
abstract type F <: Imperial end
TemperatureUnits = Union{Kelvin,K,Celsius,C,F}
default_unit(::Type{T}) where {T<:TemperatureUnits} = Kelvin
convert_unit(::Type{Kelvin},::Type{K},val) = val
convert_unit(::Type{K},::Type{Kelvin},val) = val
convert_unit(::Type{Kelvin},::Type{Celsius},val) = val + 273.15
convert_unit(::Type{Celsius},::Type{Kelvin},val) = val - 273.15
convert_unit(::Type{Kelvin},::Type{C},val) = val + 273.15
convert_unit(::Type{C},::Type{Kelvin},val) = val - 273.15
convert_unit(::Type{Kelvin},::Type{F},val) = (val-32)*5/9 + 273.15
convert_unit(::Type{F},::Type{Kelvin},val) = (val-273.15)*9/5 + 32


# pressure
abstract type Pascals <: SI end
abstract type Pa <: SIOtherUnits end
abstract type KPa <: SIOtherUnits end
abstract type atm <: OtherUnits end
abstract type psi <: OtherUnits end
PressureUnits = Union{Pascals,Pa,KPa,atm,psi}
default_unit(::Type{T}) where {T<:PressureUnits} = Pascals
convert_unit(::Type{Pascals},::Type{Pa},val) = val
convert_unit(::Type{Pa},::Type{Pascals},val) = val
convert_unit(::Type{Pascals},::Type{KPa},val) = val*1e3
convert_unit(::Type{KPa},::Type{Pascals},val) = val*1e-3
convert_unit(::Type{Pascals},::Type{atm},val) = val*101325
convert_unit(::Type{atm},::Type{Pascals},val) = val/101325
convert_unit(::Type{Pascals},::Type{psi},val) = val*6894.76
convert_unit(::Type{psi},::Type{Pascals},val) = val/6894.76


# density
abstract type KGPerCuM <: SI end
DensityUnits = Union{KGPerCuM}
default_unit(::Type{T}) where {T<:DensityUnits} = KGPerCuM

# energies
abstract type JPerKG <: SI end
SpecificEnergyUnits = Union{JPerKG}
default_unit(::Type{SpecificEnergyUnits}) = JPerKG

# enthalpy
EnthalpyUnits = SpecificEnergyUnits
default_unit(::Type{T}) where {T<:EnthalpyUnits} = default_unit(SpecificEnergyUnits)

# internal energy
InternalEnergyUnits = SpecificEnergyUnits
default_unit(::Type{T}) where {T<:InternalEnergyUnits} = default_unit(SpecificEnergyUnits)

# velocities
abstract type MPerSec <: SI end
VelocityUnits = Union{MPerSec}
default_unit(::Type{VelocityUnits}) = MPerSec

# speed of sound
SoundSpeedUnits = VelocityUnits
default_unit(::Type{T}) where {T<:SoundSpeedUnits} = default_unit(VelocityUnits)

# duct area
abstract type SqM <: SI end
abstract type SqCM <: SIOtherUnits end
AreaUnits = Union{SqM,SqCM}
default_unit(::Type{T}) where {T<:AreaUnits} = SqM
convert_unit(::Type{SqM},::Type{SqCM},val) = val*1e-4
convert_unit(::Type{SqCM},::Type{SqM},val) = val*1e4


# Mach number
MachNumberUnits = Dimensionless
default_unit(::Type{T}) where {T<:MachNumberUnits} = Dimensionless

# specific heats and gas constant
abstract type JPerKGK <: SI end
GasConstantUnits = Union{JPerKGK}
default_unit(::Type{T}) where {T<:GasConstantUnits} = JPerKGK

# entropy
EntropyUnits = Union{JPerKGK}
default_unit(::Type{T}) where {T<:EntropyUnits} = JPerKGK


###### THERMODYNAMIC PROCESSES #######

abstract type ThermodynamicProcess end

abstract type Isentropic <: ThermodynamicProcess end
abstract type NormalShock <: ThermodynamicProcess end

###### THERMODYNAMIC QUANTITIES #######

abstract type ThermodynamicQuantity{U<:ThermodynamicUnits} end

units(::ThermodynamicQuantity{U}) where {U <: ThermodynamicUnits} = U
value(s::ThermodynamicQuantity) = s.val
value(s::ThermodynamicQuantity{U},units::Type{T}) where {U <: ThermodynamicUnits, T <:ThermodynamicUnits} = convert_unit(units,U,s.val)
name(s::ThermodynamicQuantity) = s.name

function Base.show(io::IO, m::MIME"text/plain", s::ThermodynamicQuantity{U}) where {U}
    units = U == Dimensionless ? "" : U
    print(io,"$(s.name) = $(round(s.val,sigdigits=6)) $units")
end

for op in (:(+),:(-))
    @eval $op(s1::ThermodynamicQuantity{U1},s2::ThermodynamicQuantity{U1}) where {U1 <: ThermodynamicUnits} = $op(s1.val,s2.val)
    @eval $op(s::ThermodynamicQuantity,C::Real) = $op(s.val,C)
    @eval $op(C::Real,s::ThermodynamicQuantity) = $op(C,s.val)
end

for op in (:(*),:(/))
    @eval $op(s1::ThermodynamicQuantity,s2::ThermodynamicQuantity) = $op(s1.val,s2.val)
    @eval $op(s::ThermodynamicQuantity,C::Real) = $op(s.val,C)
    @eval $op(C::Real,s::ThermodynamicQuantity) = $op(C,s.val)
end

for op in (:(^),)
    @eval $op(s::ThermodynamicQuantity,C::Real) = $op(s.val,C)
end

abstract type ThermodynamicProperty{U<:ThermodynamicUnits} <: ThermodynamicQuantity{U} end

struct SpecificHeatRatio{Dimensionless} <: ThermodynamicProperty{Dimensionless}
    val :: Float64
    name :: String
end
SpecificHeatRatio(val) = SpecificHeatRatio{Dimensionless}(val,"SpecificHeatRatio")

###### THERMODYNAMIC PROPERTIES #######

for qty in ("SpecificHeatPressure","SpecificHeatVolume","GasConstant")
    qtysym = Symbol(qty)

    @eval struct $qtysym{GasConstantUnits} <: ThermodynamicProperty{GasConstantUnits}
        val :: Float64
        name :: String
    end

    @eval $qtysym(val::Real;units::Type{T}=default_unit(GasConstantUnits)) where {T<:GasConstantUnits} = $qtysym{units}(val,$qty)
end

# default values
const DefaultGasConstant = GasConstant(287)
const DefaultSpecificHeatRatio = SpecificHeatRatio(1.4)

SpecificHeatPressure(;γ::SpecificHeatRatio=DefaultSpecificHeatRatio,R::GasConstant=DefaultGasConstant) = SpecificHeatPressure(γ*R/(γ-1))
SpecificHeatVolume(;γ::SpecificHeatRatio=DefaultSpecificHeatRatio,R::GasConstant=DefaultGasConstant) = SpecificHeatVolume(R/(γ-1))

######## GAS DEFINITIONS #######

abstract type Gas end

struct PerfectGas <: Gas
    γ :: SpecificHeatRatio
    R :: GasConstant
    cp :: SpecificHeatPressure
    cv :: SpecificHeatVolume
end

function PerfectGas(;γ::SpecificHeatRatio=DefaultSpecificHeatRatio,R::GasConstant=DefaultGasConstant)
    return PerfectGas(γ,R,SpecificHeatPressure(γ=γ,R=R),SpecificHeatVolume(γ=γ,R=R))
end

SpecificHeatRatio(g::PerfectGas) = g.γ
GasConstant(g::PerfectGas) = g.R
SpecificHeatPressure(g::PerfectGas) = g.cp
SpecificHeatVolume(g::PerfectGas) = g.cv

function Base.show(io::IO, m::MIME"text/plain", g::PerfectGas)
    println(io,"Perfect gas with")
    println(io,"   Specific heat ratio = $(round(value(SpecificHeatRatio(g)),sigdigits=6))")
    println(io,"   Gas constant = $(round(value(GasConstant(g)),sigdigits=6))")
    println(io,"   cp = $(round(value(SpecificHeatPressure(g)),sigdigits=6))")
    println(io,"   cv = $(round(value(SpecificHeatVolume(g)),sigdigits=6))")

end

const DefaultPerfectGas = PerfectGas()
const Air = PerfectGas()
const He = PerfectGas(γ=SpecificHeatRatio(5/3),R=GasConstant(2077))
const O2 = PerfectGas(γ=SpecificHeatRatio(1.4),R=GasConstant(260))
const CO2 = PerfectGas(γ=SpecificHeatRatio(1.3),R=GasConstant(189))
const H2 = PerfectGas(γ=SpecificHeatRatio(1.40),R=GasConstant(4124))
const N2 = PerfectGas(γ=SpecificHeatRatio(1.40),R=GasConstant(297))

######## THERMODYNAMIC STATES #######

abstract type ThermodynamicStateVar{U<:ThermodynamicUnits} <: ThermodynamicQuantity{U} end


# thermodynamic state quantities with stagnation values
for qty in ("Pressure","Temperature","Density","Enthalpy","InternalEnergy","SoundSpeed")

    qtysym = Symbol(qty)
    stagqty = "Stagnation"*qty
    stagqtysym = Symbol(stagqty)
    unitsym = Symbol(qty,"Units")
    @eval struct $qtysym{U <: $unitsym} <: ThermodynamicStateVar{U}
        val :: Float64
        name :: String
    end

    @eval struct $stagqtysym{U <: $unitsym} <: ThermodynamicStateVar{U}
        val :: Float64
        name :: String
    end

    @eval $qtysym(val::Real;units::Type{T}=default_unit($unitsym)) where {T<:$unitsym} =
            $qtysym{default_unit($unitsym)}(convert_unit(default_unit($unitsym),units,val),$qty)

    @eval $stagqtysym(val::Real;units::Type{T}=default_unit($unitsym)) where {T<:$unitsym} =
            $stagqtysym{default_unit($unitsym)}(convert_unit(default_unit($unitsym),units,val),$stagqty)


end

# other thermodynamic state quantities
for qty in ("Area","MachNumber","Entropy")

    qtysym = Symbol(qty)
    unitsym = Symbol(qty,"Units")
    @eval struct $qtysym{U <: $unitsym} <: ThermodynamicStateVar{U}
        val :: Float64
        name :: String
    end

    @eval $qtysym(val::Real;units::Type{T}=default_unit($unitsym)) where {T<:$unitsym} =
            $qtysym{default_unit($unitsym)}(convert_unit(default_unit($unitsym),units,val),$qty)
end

# ratios of thermodynamic quantities
for qtybase in ("Area","Pressure","Temperature","Density","MachNumber","StagnationPressure")

    qty = qtybase*"Ratio"
    qtysym = Symbol(qty)
    unitsym = Symbol(qty,"Units")

    @eval $unitsym = Dimensionless

    @eval default_unit(::Type{T}) where {T<:$unitsym} = Dimensionless

    @eval struct $qtysym{U <: $unitsym} <: ThermodynamicStateVar{U}
        val :: Float64
        name :: String
    end

    @eval $qtysym(val::Real;units::Type{T}=default_unit($unitsym)) where {T<:$unitsym} =
            $qtysym{default_unit($unitsym)}(convert_unit(default_unit($unitsym),units,val),$qty)
end

#=
Thermodynamic relationships
=#

####### BASIC PERFECT GAS EQUATIONS OF STATE #######

Temperature(ρ::Density,p::Pressure;gas::PerfectGas=DefaultPerfectGas) = Temperature(p/(ρ*GasConstant(gas)))
Density(p::Pressure,T::Temperature;gas::PerfectGas=DefaultPerfectGas) = Density(p/(GasConstant(gas)*T))
Pressure(ρ::Density,T::Temperature;gas::PerfectGas=DefaultPerfectGas) = Pressure(ρ*GasConstant(gas)*T)

"""
    SoundSpeed(T::Temperature[,γ::SpecificHeatRatio=1.4][,R::GasConstant=287])

Compute the speed of sound of a perfect gas, based on the absolute temperature `T`. The ratio of specific heats `γ` and `R`
can also be supplied, but default to the values for air at standard conditions
(1.4 and 287 J/kg.K, respectively).
"""
SoundSpeed(T::Temperature;gas::PerfectGas=DefaultPerfectGas) = SoundSpeed(sqrt(SpecificHeatRatio(gas)*GasConstant(gas)*T))

Enthalpy(T::Temperature;gas::PerfectGas=DefaultPerfectGas) = Enthalpy(SpecificHeatPressure(gas)*T)
StagnationEnthalpy(T0::StagnationTemperature;gas::PerfectGas=DefaultPerfectGas) = StagnationEnthalpy(SpecificHeatPressure(gas)*T0)

InternalEnergy(T::Temperature;gas::PerfectGas=DefaultPerfectGas) = InternalEnergy(SpecificHeatVolume(gas)*T)
StagnationInternalEnergy(T0::StagnationTemperature;gas::PerfectGas=DefaultPerfectGas) = StagnationInternalEnergy(SpecificHeatVolume(gas)*T0)


function T0OverT(M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    return TemperatureRatio(1+0.5*(γ-1)*M^2)
end

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

function DensityRatio(Tratio::TemperatureRatio,::Type{Isentropic};gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    DensityRatio(Tratio^(1/(γ-1)))
end

# Temperature relations
function StagnationTemperature(T::Temperature,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    StagnationTemperature(T*T0OverT(M,gas=gas))
end

Temperature(T0::StagnationTemperature,M::MachNumber;gas::PerfectGas=DefaultPerfectGas) = Temperature(T0/T0OverT(M,gas=gas))

function MachNumber(T_over_T0::TemperatureRatio;gas::PerfectGas=DefaultPerfectGas)
    value(T_over_T0) <= 1.0 || error("T/T0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    M2 = ((1.0/T_over_T0)-1)*2/(γ-1)
    MachNumber(sqrt(M2))
end

MachNumber(T::Temperature,T0::StagnationTemperature;gas::PerfectGas=DefaultPerfectGas) = MachNumber(TemperatureRatio(T/T0),gas=gas)

# Pressure relations
function StagnationPressure(p::Pressure,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    StagnationPressure(p*PressureRatio(T0OverT(M,gas=gas),Isentropic))
end
function Pressure(p0::StagnationPressure,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    Pressure(p0/PressureRatio(T0OverT(M,gas=gas),Isentropic))
end

function MachNumber(p_over_p0::PressureRatio;gas::PerfectGas=DefaultPerfectGas)
    value(p_over_p0) <= 1.0 || error("p/p0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    MachNumber(TemperatureRatio(p_over_p0,Isentropic,gas=gas))
end

MachNumber(p::Pressure,p0::StagnationPressure;gas::PerfectGas=DefaultPerfectGas) = MachNumber(PressureRatio(p/p0),gas=gas)

# Density relations
function StagnationDensity(ρ::Density,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    StagnationDensity(ρ*DensityRatio(T0OverT(M,gas=gas),Isentropic))
end
function Density(ρ0::StagnationDensity,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    Density(ρ0/DensityRatio(T0OverT(M,gas=gas),Isentropic))
end

function MachNumber(ρ_over_ρ0::DensityRatio;gas::PerfectGas=DefaultPerfectGas)
    value(ρ_over_ρ0) <= 1.0 || error("ρ/ρ0 must be 1 or smaller")
    γ = SpecificHeatRatio(gas)
    MachNumber(TemperatureRatio(ρ_over_ρ0,Isentropic,gas=gas))
end

MachNumber(ρ::Density,ρ0::StagnationDensity;gas::PerfectGas=DefaultPerfectGas) = MachNumber(DensityRatio(ρ/ρ0),gas=gas)

# Area-Mach number relation
function AOverAStar(M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    γ = SpecificHeatRatio(gas)
    return AreaRatio(1/M*(2/(γ+1)*(T0OverT(M,gas=gas)))^(0.5*(γ+1)/(γ-1)))
end

function AStar(A::Area,M::MachNumber;gas::PerfectGas=DefaultPerfectGas)
    return Area(A/AOverAStar(M,gas=gas))
end

function MachNumber(A_over_Astar::AreaRatio;gas::PerfectGas=DefaultPerfectGas)
    value(A_over_Astar) >= 1.0 || error("A/A* must be 1 or larger")
    Msub = find_zero(x -> AOverAStar(MachNumber(x))-A_over_Astar,(0,1),order=16)
    Msup = find_zero(x -> AOverAStar(MachNumber(x))-A_over_Astar,(1,Inf),order=16)
    return MachNumber(Msub), MachNumber(Msup)
end

function MachNumber(M1::MachNumber,A1::Area,A2::Area;gas::PerfectGas=DefaultPerfectGas)
    M2sub, M2sup = MachNumber(AreaRatio(A2/AStar(A1,M1,gas=gas)),gas=gas)
end

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
    p01_unit = StagnationPressure(Pressure(1),M1,gas=gas)
    p02_unit = StagnationPressure(Pressure(1),M2,gas=gas)
    return StagnationPressureRatio(pratio*p02_unit/p01_unit)
end

function Entropy(s1::Entropy,M1::MachNumber,::Type{NormalShock};gas::PerfectGas=DefaultPerfectGas)
    value(M1) >= 1.0 || error("M1 must be at least 1")
    γ = SpecificHeatRatio(gas)
    cv = SpecificHeatVolume(gas)
    ds = cv*(log(1+2*γ/(γ+1)*(M1^2-1)) - γ*log((γ+1)*M1^2/(2+(γ-1)*M1^2)))
    return Entropy(s1+ds)

end

nothing
