include("header.jl")

#=
# Basic notes on complex variables
This notebook reviews some aspects of complex variables and demonstrates how to use
them in Julia.
=#

# ### Set up the module
using MAE150A
#-
using Plots

#=
Potential flows make common use of **complex variables**. For example, rather than
write the position of a point in Cartesian coordinates as $(x,y)$, it is typical
to write it instead as

$$ z = x + \mathrm{i} y $$

where $\mathrm{i} = \sqrt{-1}$ is the imaginary unit. It is also common for us to use
cylindrical (or polar) coordinates. These also have a complex form. If a point is at
$(r,\theta)$, then we write this as

$$ z = r \mathrm{e}^{\mathrm{i}\theta} $$

Remember Euler's formula, $\mathrm{e}^{\mathrm{i}\theta} = \cos\theta + \mathrm{i}\sin\theta$.
This ensures that $x = r \cos\theta$ and $y = r \sin\theta$, as we should expect.
=#

#=
#### Complex variables in Julia
Julia allows us to write complex variables very easily. The imaginary unit is `im`.
So, for example, the point $(1.0,2.0)$ can be written in complex form as
=#

z = 1.0+im*2.0


# (The reason we used 1.0 instead of 1 was to ensure that this is treated as a floating-point
# number, rather than an integer.)
#
# We can always recover the real or imaginary part of a variable with `real` and `imag`:

real(z)
#-
imag(z)

# We can also set up an array of complex numbers. For example, here is a set of points at
# $y = 0.1$ and between $x = -1$ and $x = 1$.
z = (-1:0.1:1) .+ 0.1im;

#=
The `.+` signifies to add `0.1im` (i.e., $0.1\mathrm{i}$) to every element in the list
`-1:0.1:1`.

If we want just the $x$ components of this array of points, then use `real`, but you will
need to use a `.` to apply it to every element of the array:
=#

real.(z)

# Plot these points just to see what we've created:
scatter(real.(z),imag.(z),ratio=1,xlim=(-2,2),ylim=(-1,1),legend=false,xlabel=L"x",ylabel=L"y")


#=
Suppose we want points on a circle. That's easy to do in the polar form,
$z = a \mathrm{e}^{\mathrm{i}\theta}$, where $a$ is the radius of the circle, and $\theta$
varies from 0 to $2\pi$.
Set the radius of the circle
=#
a = 0.5

# Set a range of angles between 0 and 2π
θ = range(0,2π,length=50)

# And finally, create the coordinates of all of the points on the circle:
z = a*exp.(im*θ)

# The last element of z is identical to the first one, so we can get rid of the redundant one
# with the following:
pop!(z);

# Let us plot these points to check. The `ratio=1` argument makes sure that we get equal
# dimensions of the axes. So a circle looks like a circle rather than a squashed circle
scatter(real.(z),imag.(z),ratio=1,xlim=(-2,2),ylim=(-2,2),legend=false,xlabel=L"x",ylabel=L"y")

#=
### Other operations
The **complex conjugate** is an important operation, denoted by $()^{*}$. For example,
if $z = x + \mathrm{i}y$, then the conjugate simply switches the sign of the imaginary part.

$$ z^{*} = x - \mathrm{i}y $$

It does not look that important, but it is useful for many operations.

For example, the magnitude squared (i.e., the 'length' squared) of a complex variable
can be obtained with the help of the conjugate, $|z|^2 = z z^* = x^2 + y^2$. Check this:

$$ zz^* = (x+\mathrm{i}y)(x-\mathrm{i}y) = x^2 + \mathrm{i}yx - \mathrm{i}xy - \mathrm{i}^2 y^2$$.

But since $\mathrm{i}^2 = -1$, then this is simply $x^2+y^2$.
=#

#=
#### Complex conjugate of the polar form
The complex conjugate of the polar form is also easy:

$$ z^* = r \mathrm{e}^{-\mathrm{i}\theta}$$.

Then, using the polar form, the magnitude squared is

$$ |z|^2 = z z^* = r \mathrm{e}^{\mathrm{i}\theta} r \mathrm{e}^{-\mathrm{i}\theta} = r^2$$

As expected, $r^2 = x^2 + y^2$.
=#

# In Julia, we use `conj` to get the complex conjugate and `abs` to get the magnitude
z = 1.0+im*2.0
#-
conj(z)
#-
abs(z)

#=
### The rotation operation
The complex number $\mathrm{e}^{\mathrm{i}\alpha}$ has special importance: By multiplying
it with another complex number, it **rotates** that number by angle $\alpha$ in the
counter-clockwise direction about the origin. That can be particularly useful for moving
points around.

Let's try this:
=#
α = π/2              ## angle of counter-clockwise rotation
z = 1.0+2.0*im       ## original point
zrot = exp(im*α)*z;  ## new point, rotated

#=
Plot these points, using some lines from the origin to show that the rotation
preserves the length.
=#
scatter(real.([z,zrot]),imag.([z,zrot]),xlim=(-3,3),ylim=(-3,3),legend=true,ratio=1)
plot!([0,real(z)],[0,imag(z)],label="Ray to original pt")
plot!([0,real(zrot)],[0,imag(zrot)],label="Ray to rotated pt")

#=
### Other uses of complex variables in potential flow
There are many other uses of complex variables in potential flow. We don't need to
cover them all here. However, it is helpful to know that we can write a **velocity vector**
in complex form in terms of its Cartesian components:

$$ u + \mathrm{i} v $$

We can also write velocity in its polar components:

$$ u_r + \mathrm{i} u_\theta $$

We can go from one form to the other by using the **local** rotation operator,
$\mathrm{e}^{-\mathrm{i}\theta}$.

$$ u_r + \mathrm{i} u_\theta = (u + \mathrm{i} v) \mathrm{e}^{-\mathrm{i}\theta}$$"
=#

# Try this rotation:
u = 0.0
v = 1.0
θ = π/2
(u+im*v)*exp(-im*θ) # ur + i*uθ

# The real and imaginary parts of this are $u_r$ and $u_\theta$
