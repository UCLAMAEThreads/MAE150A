include("header.jl")

#=
# Basic Tools for Quasi-1D Steady Compressible Flow
This notebook demonstrates the basic syntax for some tools for computing
quasi-1d steady compressible flow.
=#

# ### Set up the module
using MAE150A
#-
using Plots

#=
### Setting basic properties and states
We can set thermodynamic properties and states in a straightforward manner.
However, it is important to remember that we have to explicitly define the
type of property or state we are setting. Examples below will show how this works.
=#

#=
#### Units
When we set a thermodynamic quantity, we also need to set its **units**. For every
quantity, we have a few units to choose from. For example:
=#
PressureUnits

#=
Among these, one is designated the default unit (the SI unit of that quantity).
All quantities are converted to this unit for calculation purposes. For example,
for pressure:
=#
default_unit(PressureUnits)

#=
or enthalpy:
=#
default_unit(EnthalpyUnits)

#=
If we do not specify the unit, it is **automatically set to the default unit**:
=#
Pressure(10000)

#=
If we set the pressure in another unit, it will still convert it to the default
unit for us. This ensures that all calculations are carried out in standardized fashion.
=#
p = Pressure(1,units=atm)

#=
However, we can always report the quantity in some desired units with the `value`
function:
=#
value(p,atm)
#-
value(p,psi)
#-
value(p,KPa)

#=
#### Other thermodynamic quantities
We can set most any other thermodynamic quantity in similar fashion:
=#
T = Temperature(20,units=Celsius)
#-
T0 = StagnationTemperature(20)
#-
MachNumber(2.0)
#-
Enthalpy(50)
#-
Entropy(10)
#-
Area(50,units=SqCM)
#-
Length(5)
#=
and others...
=#

#=
#### Gas properties
We can set the properties of the gas that we are analyzing. (Note: It is
assumed that the gas is perfect.)
=#
SpecificHeatRatio(1.3)
#-
GasConstant(320)

#=
and we can define a gas with these values:
=#
gas = PerfectGas(γ=SpecificHeatRatio(1.3),R=GasConstant(320))

#=
We have **pre-defined gases** (at standard conditions), as well, for convenience:
=#
Air
#-
He
#-
O2
#-
CO2
#-
H2
#-
N2
#=
#### Equations of state
We can apply the equation of state for a perfect gas to determine other quantities.
For example, suppose we have carbon dioxide at 1.2 kg/m^3 and 80 KPa. What is the temperature?
=#
T = Temperature(Density(1.2),Pressure(80,units=KPa),gas=CO2)

#=
You can switch the order of the arguments and it will still work:
=#
T = Temperature(Pressure(80,units=KPa),Density(1.2),gas=CO2)

#=
Then we can calculate the enthalpy, for example:
=#
Enthalpy(T,gas=CO2)

#=
What is the speed of sound of air at 20 degrees Celsius? Let's find out:
=#
SoundSpeed(Temperature(20,units=C),gas=Air)

#=
How about oxygen?
=#
SoundSpeed(Temperature(20,units=C),gas=O2)

#=
**Note: the default gas is air. So if you do not put the `gas=` argument in,
it will assume air at standard conditions.**
=#
