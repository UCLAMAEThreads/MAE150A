include("header.jl")

#=
# Basic Tools for Quasi-1D Steady Compressible Flow
This notebook demonstrates the basic syntax for some tools for computing
quasi-1d steady compressible flow.
=#

# ### Set up the module
using MAE150A
#-
using Gasdynamics1D

#=
### Setting basic properties and states
We can set thermodynamic properties and states in a straightforward manner.
However, it is important to remember that we have to explicitly define the
type of property or state we are setting. Examples below will show how this works.
=#

#=
Let's say we wish to set the pressure to 10000 Pa. Pascals are the default units
of pressure, as can be verified by using the `default_unit` function:
=#
default_unit(Pressure)

#=
So if we do not specify the unit, it is **automatically set to the default unit**:
=#
Pressure(10000)


#=
We can set a quantity with another unit using the syntax u"[unit]". For example,
if we set pressure to 1 atm, it will still convert it to the default
unit.
=#
p = Pressure(1u"atm")

#=
However, we can always report the quantity in some desired units with the `value`
function:
=#
value(p,u"atm")
#-
value(p,u"psi")
#-
value(p,u"kPa")

#=
#### Other thermodynamic quantities
We can set most any other thermodynamic quantity in similar fashion:
=#
T = Temperature(20u"°C")
#-
T0 = StagnationTemperature(20)
#-
MachNumber(2.0)
#-
Enthalpy(50)
#-
Entropy(10)
#-
Area(50u"cm^2")
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
#-
Ar
#=
#### Equations of state
We can apply the equation of state for a perfect gas to determine other quantities.
For example, suppose we have carbon dioxide at 1.2 kg/m^3 and 80 kPa. What is the temperature?
=#
T = Temperature(Density(1.2),Pressure(80u"kPa"),gas=CO2)

#=
You can switch the order of the arguments and it will still work:
=#
T = Temperature(Pressure(80u"kPa"),Density(1.2),gas=CO2)

#=
Then we can calculate the enthalpy, for example:
=#
Enthalpy(T,gas=CO2)

#=
What is the speed of sound of air at 20 degrees Celsius? Let's find out:
=#
SoundSpeed(Temperature(20u"°C"),gas=Air)

#=
How about oxygen?
=#
SoundSpeed(Temperature(20u"°C"),gas=O2)

#=
**Note: the default gas is air. So if you do not put the `gas=` argument in,
it will assume air at standard conditions.**
=#
