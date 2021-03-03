#=
Quasi-1d gas dynamics routines
=#

import Base:+,*,-,/,^,>,<,>=,<=,==,isapprox

using Roots

export default_unit,ThermodynamicUnits,SI,SIOtherUnits,Imperial,OtherUnits,Dimensionless
export convert_unit
export Kelvin, K, Celsius, C, F, TemperatureUnits
export Pascals, Pa, KPa, atm, psi, PressureUnits
export KGPerCuM, DensityUnits
export JPerKG, KJPerKG, SpecificEnergyUnits, EnthalpyUnits, InternalEnergyUnits,
        HeatFluxUnits
export MPerSec, VelocityUnits, SoundSpeedUnits
export SqM, SqCM,AreaUnits
export MachNumberUnits
export JPerKGK, KJPerKGK, GasConstantUnits
export EntropyUnits
export KGPerSec, MassFlowRateUnits
export Meters, CM, LengthUnits, DiameterUnits
export FrictionFactorUnits, FLOverDUnits

export ThermodynamicProcess, Isentropic, NormalShock, FannoFlow, RayleighFlow
export ThermodynamicQuantity
export units, value, name
export ThermodynamicProperty
export SpecificHeatRatio
export SpecificHeatPressure,SpecificHeatVolume,GasConstant
export DefaultGasConstant, DefaultSpecificHeatRatio
export Gas, PerfectGas, DefaultPerfectGas
export Air, He, O2,CO2,H2,N2
export ThermodynamicStateVar, Pressure, Temperature, Density, Enthalpy, InternalEnergy, SoundSpeed
export Area, MachNumber, Entropy, MassFlowRate,Velocity,Length, Diameter
export SubsonicMachNumber,SupersonicMachNumber
export AreaRatio, FrictionFactor, FLOverD, HeatFlux
export StagnationPressure, StagnationDensity, StagnationEnthalpy,StagnationTemperature
export StagnationPressureRatio, PressureRatio, DensityRatio, TemperatureRatio
export T0OverT,P0OverP,ρ0Overρ,AOverAStar,AStar
export FLStarOverD,POverPStar,ρOverρStar,TOverTStar,P0OverP0Star
export HeatFlux,T0OverT0Star,VOverVStar



#=
Note: to add a new quantity, you define its units here, calling its units, for
example, MyQtyUnits. Then, add MyQty below to the list of strings in
Thermodynamic States section below.
=#

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
abstract type KJPerKG <: SIOtherUnits end
SpecificEnergyUnits = Union{JPerKG,KJPerKG}
default_unit(::Type{SpecificEnergyUnits}) = JPerKG
convert_unit(::Type{JPerKG},::Type{KJPerKG},val) = val*1e3
convert_unit(::Type{KJPerKG},::Type{JPerKG},val) = val*1e-3

# enthalpy
EnthalpyUnits = SpecificEnergyUnits
#default_unit(::Type{T}) where {T<:EnthalpyUnits} = default_unit(SpecificEnergyUnits)

# internal energy
InternalEnergyUnits = SpecificEnergyUnits
#default_unit(::Type{T}) where {T<:InternalEnergyUnits} = default_unit(SpecificEnergyUnits)

# heat flux
HeatFluxUnits = SpecificEnergyUnits
#default_unit(::Type{T}) where {T<:HeatFluxUnits} = default_unit(SpecificEnergyUnits)

# velocities
abstract type MPerSec <: SI end
VelocityUnits = Union{MPerSec}
default_unit(::Type{VelocityUnits}) = MPerSec



# speed of sound
SoundSpeedUnits = VelocityUnits
#default_unit(::Type{T}) where {T<:SoundSpeedUnits} = default_unit(VelocityUnits)



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
abstract type KJPerKGK <: SIOtherUnits end
GasConstantUnits = Union{JPerKGK,KJPerKGK}
default_unit(::Type{T}) where {T<:GasConstantUnits} = JPerKGK
convert_unit(::Type{JPerKGK},::Type{KJPerKGK},val) = val*1e3
convert_unit(::Type{KJPerKGK},::Type{JPerKGK},val) = val*1e-3


# entropy
EntropyUnits = GasConstantUnits
#default_unit(::Type{T}) where {T<:EntropyUnits} = JPerKGK

# mass flow rate
abstract type KGPerSec <: SI end
MassFlowRateUnits = Union{KGPerSec}
default_unit(::Type{T}) where {T<:MassFlowRateUnits} = KGPerSec

# length
abstract type Meters <: SI end
abstract type CM <: SIOtherUnits end
LengthUnits = Union{Meters,CM}
default_unit(::Type{LengthUnits}) = Meters
convert_unit(::Type{Meters},::Type{CM},val) = val*1e-2
convert_unit(::Type{CM},::Type{Meters},val) = val*1e2


# diameter
DiameterUnits = LengthUnits
#default_unit(::Type{T}) where {T<:DiameterUnits} = default_unit(LengthUnits)

# friction factor f
FrictionFactorUnits = Dimensionless
#default_unit(::Type{T}) where {T<:FrictionFactorUnits} = Dimensionless

# f*L/D
FLOverDUnits = Dimensionless
#default_unit(::Type{T}) where {T<:FLOverDUnits} = Dimensionless


###### THERMODYNAMIC PROCESSES #######

abstract type ThermodynamicProcess end

abstract type Isentropic <: ThermodynamicProcess end
abstract type NormalShock <: ThermodynamicProcess end
abstract type FannoFlow <: ThermodynamicProcess end
abstract type RayleighFlow <: ThermodynamicProcess end


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

for op in (:(+),:(-),:(>),:(<),:(>=),:(<=),:(==))
    @eval $op(s1::ThermodynamicQuantity{U1},s2::ThermodynamicQuantity{U1}) where {U1 <: ThermodynamicUnits} = $op(s1.val,s2.val)
    @eval $op(s::ThermodynamicQuantity,C::Real) = $op(s.val,C)
    @eval $op(C::Real,s::ThermodynamicQuantity) = $op(C,s.val)
end

for op in (:(isapprox),)
    @eval $op(s1::ThermodynamicQuantity{U1},s2::ThermodynamicQuantity{U1};kwargs...) where {U1 <: ThermodynamicUnits} = $op(s1.val,s2.val;kwargs...)
    @eval $op(s::ThermodynamicQuantity,C::Real;kwargs...) = $op(s.val,C;kwargs...)
    @eval $op(C::Real,s::ThermodynamicQuantity;kwargs...) = $op(C,s.val;kwargs...)
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
const H2 = PerfectGas(γ=SpecificHeatRatio(1.405),R=GasConstant(4126))
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
for qty in ("Area","MachNumber","Entropy","MassFlowRate","Velocity",
            "Length","Diameter","FrictionFactor","FLOverD","HeatFlux")

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
for qtybase in ("Area","Pressure","Temperature","Density","MachNumber",
                "StagnationPressure","Velocity")

    qty = qtybase*"Ratio"
    qtysym = Symbol(qty)
    unitsym = Symbol(qty,"Units")

    @eval $unitsym = Dimensionless

    #@eval default_unit(::Type{T}) where {T<:$unitsym} = Dimensionless

    @eval struct $qtysym{U <: $unitsym} <: ThermodynamicStateVar{U}
        val :: Float64
        name :: String
    end

    @eval $qtysym(val::Real;units::Type{T}=default_unit($unitsym)) where {T<:$unitsym} =
            $qtysym{default_unit($unitsym)}(convert_unit(default_unit($unitsym),units,val),$qty)
end

#=
Geometric relationships
=#
Area(D::Diameter) = Area(π*D^2/4)
Diameter(A::Area) = Diameter(sqrt(4*A/π))



include("gasdynamics/thermodynamics.jl")
include("gasdynamics/isentropic.jl")
include("gasdynamics/normalshocks.jl")
include("gasdynamics/fanno.jl")
include("gasdynamics/rayleigh.jl")
include("gasdynamics/nozzle.jl")
include("gasdynamics/vectors.jl")
