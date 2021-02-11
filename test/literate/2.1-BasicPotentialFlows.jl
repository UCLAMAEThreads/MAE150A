include("header.jl")

#=
# Basic planar potential flows
In this notebook we will discuss elementary potential flows. These will be the building
blocks of more complex potential flows that come later.

Remember that, by definition, all of these flows are solutions of Laplace's equation,
both for the scalar potential

$$ \nabla^2 \phi = 0$$

and the streamfunction

$$ \nabla^2 \psi = 0$$

In some cases (singularities), however, they fail to solve these equations at some points.

Also, remember the relationships of these quantities with velocity:

$$ u = \dfrac{\partial \phi}{\partial x} = \dfrac{\partial \psi}{\partial y}$$

and

$$ v = \dfrac{\partial \phi}{\partial y} = -\dfrac{\partial \psi}{\partial x}$$

In polar coordinates,

$$ u_r = \dfrac{\partial \phi}{\partial r} = \dfrac{1}{r} \dfrac{\partial \psi}{\partial \theta}$$

and

$$ u_\theta = \dfrac{1}{r} \dfrac{\partial \phi}{\partial \theta} = -\dfrac{\partial \psi}{\partial r}$$
=#

# ### Set up the module
using MAE150A
#-
using PotentialFlow
using ViscousFlow
#-
using Plots

#=
### Set up grids and grid data for visualizing the potential flows
To see these flows, we need to evaluate them on a grid of points. As we
discussed in the previous notebook, it is most convenient to evaluate these
using complex coordinates. So here, we will set up grids of complex coordinates
for evaluating the different flow types.
=#

# First make the grid itself
## grid cell size
Δx = 0.02

## set the size of a region in which to set up the flow.
xlim = (-2,2)
ylim = (-2,2)

## make the grid
g = PhysicalGrid(xlim,ylim,Δx)

# And now set up streamfunction and velocity and their associated grid coordinates
ψ = Nodes(Dual,size(g)) # streamfunction field
u = Edges(Primal,size(g)) # velocity field
xg, yg = coordinates(ψ,g)
xu, yu, xv, yv = coordinates(u,g);

# and let's set up complex versions of this grid for evaluating the flows
zg = complexgrid(xg,yg)
zu = complexgrid(xu,yu)
zv = complexgrid(xv,yv);

#=
Hang on... what? Complex variables? Make sure to review the notes in the
[complex variables review notebook](2.1-ComplexVariablesNotes.ipynb) if you want to know
why we would do this or if you feel uncomfortable with your complex variable knowledge...

Let's look at the complex form of the grid points we just set up:
=#
zg

#=
Each of the points in this array is a unique $(x,y)$ location on the grid. The top
left element in the array is actually the bottom left point on the grid.
=#

# ## The basic building block flows

#=
### The most basic: Uniform flow
First, let us inspect a uniform flow, also called a *free stream*.

$$ u = U_\infty \cos\alpha,\quad v = U_\infty \sin \alpha, \qquad \psi = U_\infty (y \cos\alpha - x \sin\alpha), \qquad \phi = U_\infty (x\cos\alpha + y \sin\alpha) $$

We can specify this with the strength of the flow. Let us set up a uniform flow with speed
equal to 1 at an angle $\alpha = 45$ degrees ($\pi/4$ radians). We will use the complex
polar notation for this:
=#

U∞ = 1.0  ## speed
α = π/3  ## angle in radians (60 degrees)
fs = Freestreams.Freestream(U∞*exp(im*α))

#=
#### Evaluation of the flow
Let's look at the streamlines of this uniform flow. For that, we evaluate the streamfunction
on our grid (`zg`) that we set up earlier.

A few notes on the command below. We use the notation `PotentialFlow.streamfunction`
to let it know that we are using the streamfunction command in the `PotentialFlow` package.
Also, the `.=` makes sure to put the evaluated streamfunction on this grid into
an array we have already set up (`ψ`)."
=#
ψ .= PotentialFlow.streamfunction(zg,fs);

# Plot the streamfunction contours
p = plot(ψ,g,xlim=(-2,2),ylim=(-2,2),color=:black,xlabel=L"x",ylabel=L"y",title="Streamlines of a uniform flow",show=true)
add_arrows!(p,fs)

#=
As expected, the streamlines are angled at 60 degrees.

Try some different angles to see the result.
=#

#=
We can also evaluate the **velocity field** of the free stream, using the
`induce_velocity` function. It should show the same value everywhere. (Note
that the last argument is actually the time, but this is a steady flow,
so it is irrelevant.)
=#
induce_velocity(zg,fs,0.0)

#=
### Basic singularity: A source
Now, let's consider a source flow
$$ u_r = Q/(2\pi r), \quad u_\theta = 0,\qquad \psi = \dfrac{Q \theta}{2\pi}, \qquad \phi = \dfrac{Q }{2\pi} \ln r$$

Let's place one at the origin with strength $Q$ equal to 1.
=#
zs = 0.0+im*0.0  ## location of the source
Q = 1.0  ## strength of the source
s = Source.Point(zs,Q)

# Evaluate its streamfunction and plot it:
ψ .= PotentialFlow.streamfunction(zg,s);
p = plot(ψ,g,xlim=(-2,2),ylim=(-2,2),color=:black,xlabel=L"x",ylabel=L"y",title="Streamlines of a source",show=true)
add_arrows!(p,s)

#=
This looks as expected, but a little strange along the $-x$ axis. Remember, the
streamfunction is multi-valued. This dark line is the **branch cut** of the streamfunction,
where it jumps from one value to another across that line.
=#

#=
**SIDE NOTE:** We cannot avoid the branch cut, but we can move it to a different
ray by using a rotation operator. To move it to some specified angle $\tau$, we
rotate all of the evaluation points from $\tau$ to $-\pi$.
=#
τ = π/4 # angle at which we prefer the branch cut
rot = exp(-im*(π+τ))  # rotation operator, which moves rotates from $\tau$ to $-\pi$.
ψ .= PotentialFlow.streamfunction(zg.*rot,s);
plot(ψ,g,color=:black,xlim=(-2,2),ylim=(-2,2),xlabel=L"x",ylabel=L"y",title="Streamlines of a source")
#src add_arrows!(p2,s)


#=
Let's evaluate the velocity field of the source. We will be a bit complicated about
this, mostly for plotting purposes. Remember from the [field plotting notebook](1.2-PlottingFlowFields),
the u and v components are stored at different places on a staggered grid. We therefore
evaluate each component on a different set of points:
=#
u.u .= real.(induce_velocity(zu,s,0.0))
u.v .= imag.(induce_velocity(zv,s,0.0));
plot(
plot(u.u,g,levels=range(-1,1,length=31),clim=(-1,1),xlabel=L"x",ylabel=L"y",title="u component"),
plot(u.v,g,levels=range(-1,1,length=31),clim=(-1,1),xlabel=L"x",ylabel=L"y",title="v component")
)

#=
### Another basic singularity: A point vortex
Now, let's consider a point vortex,

$$ u_r = 0, \quad u_\theta = \dfrac{\Gamma}{2\pi r},\qquad \psi = -\dfrac{\Gamma}{2\pi} \ln r, \qquad \phi = \dfrac{\Gamma \theta}{2\pi} $$

Let's place one at the origin with strength $\Gamma$ equal to 1.
=#
zv = 0.0+im*0.0  # location of the vortex
Γ = 1.0  # strength of the vortex
v = Vortex.Point(zv,Γ)

# Evaluate its streamfunction and plot it:
ψ .= PotentialFlow.streamfunction(zg,v);
p = plot(ψ,g,color=:black,xlabel=L"x",ylabel=L"y",title="Streamlines of a vortex",show=true)
add_arrows!(p,v)

#=
### Another singularity: a dipole (or doublet)
A doublet is also a singularity, but higher order than the source and vortex. It has the form

$$ u_r = -\dfrac{D}{\pi r^2} \cos(\theta-\alpha),\quad u_\theta = -\dfrac{D}{\pi r^2} \sin(\theta-\alpha), \qquad \psi = -\dfrac{D}{\pi r} \sin(\theta-\alpha), \qquad \phi = \dfrac{D}{\pi r} \cos(\theta-\alpha) $$

$D$ is the strength and $\alpha$ the angle. Notice that it has both a $u_r$ and a $u_\theta$
component of velocity, and has a different dependence on $r$. It decays faster with distance
away from the center.
=#

# Let's create a doublet of strength 1 at angle $\pi/4$:
zd = 0.0+im*0.0
D = 1.0
α = π/4
d = Doublets.Doublet(zd,D*exp(im*α))


ψ .= PotentialFlow.streamfunction(zg,d);
p = plot(ψ,g,xlim=(-2,2),ylim=(-2,2),levels=range(-1,1,length=15),color=:black,xlabel=L"x",ylabel=L"y",title="Streamlines of a doublet",show=true)
add_arrows!(p,d)


#=
### A corner flow
A corner flow is a flow in a sector of interior angle $\nu\pi$. It is generated by the
streamfunction

$$ \psi(r,\theta) = \nu \sigma r^{1/\nu} \cos((\theta-\alpha)/\nu)$$

The strength is $\sigma$ and the rotation angle of the corner is $\alpha$.
=#

# Create a corner flow of strength 1 and interior angle $\pi/3$ (60 degrees).
ν = 0.5
σ = 1.0
α = π/3
c = Corner(σ,ν,α)

ψ .= PotentialFlow.streamfunction(zg,c);
p = plot(ψ,g,xlim=(-2,2),ylim=(-2,2),levels=range(-3,3,length=31),color=:black,xlabel=L"x",ylabel=L"y",title="Streamlines of a corner",show=true)
add_arrows!(p,c)

#=
Note the straight streamlines that cross at the origin. There is a stagnation point there.
=#

#=
## Combinations of potential flows
We can easily make combinations of potential flows, simply by adding them.

Let's try a combination of a uniform flow at 0 degrees and a source at the origin of strength 2:
=#

# The uniform flow
U∞ = 1.0  ## speed
α = 0.0 ## angle in radians
f = Freestreams.Freestream(U∞*exp(im*α))

# The source
zs = 0.0+im*0.0  ## location of the source
Q = 2.0  ## strength of the source
s = Source.Point(zs,Q)

#=
We will add these, but take some care to use our rotation trick on the branch cut
of the source, so that it is along the $+x$ axis.
=#
τ = 0.0 ## angle at which we prefer the branch cut
rot = exp(-im*(π+τ))
ψ .= PotentialFlow.streamfunction(zg,f) .+ PotentialFlow.streamfunction(zg.*rot,s);
p = plot(ψ,g,color=:black,xlim=(-2,2),ylim=(-2,2),xlabel=L"x",ylabel=L"y",title="Streamlines of a uniform flow + source",show=true)
add_arrows!(p,(f,s))

#=
There is a **stagnation point** in this flow somewhere to the left of the source,
and this indicates that **stagnation streamlines** cross there.

To include the stagnation streamline to the plot, we can find the value of streamfunction
at the stagnation point. This is the value it will have on the entire stagnation streamlines.

In this example, we can find the location of the stagnation point by hand, using the
polar velocity components from the flow:

$$ u_r = U_\infty \cos\theta + \dfrac{Q}{2\pi r}, \quad u_\theta = 0 $$

The stagnation point is at $\theta = \pi$ and $r = a = Q/(2\pi U_\infty)$. In Cartesian
coordinates, this is $x = -Q/(2\pi U_\infty)$ and $y = 0$.
=#

#=
To evaluate the streamfunction at a specific point, we first turn the grid of
streamfunction values into a function of $x$ and $y$.
=#
ψfield = interpolatable_field(xg,yg,ψ);

# Now evaluate $\psi$ at the stagnation point:
xstag = -Q/(2π*U∞)
ystag = 0
ψstag = ψfield(xstag,ystag)

# Now add a streamline with this value of $\psi$ to the plot:
p = plot(ψ,g,color=:black,xlim=(-2,2),ylim=(-2,2),xlabel=L"x",ylabel=L"y",title="Streamlines of a uniform flow + source",show=true)
plot!(p,ψ,g,levels=[ψstag],linewidth=2)
scatter!(p,[xstag],[ystag],label="stagnation point")
add_arrows!(p,(f,s))

#=
### Another combination: two vortices
The tools provided in `PotentialFlow` do not require us to explicitly add the
flows; it does it for us. This is especially useful with a lot of vortices. Let us
try two point vortices, and place them at $(-1,0)$ and $(1,0)$, and give them strengths
$1$ and $-1$:"
=#

# First create empty arrays
zvort = ComplexF64[]
Γvort = Float64[]

# Add the first vortex to the array
push!(zvort,-1.0+0im)
push!(Γvort,1.0)

# add the second vortex to the array
push!(zvort,1.0+0im)
push!(Γvort,-1.0)

#=
Now make the list of vortices. Note the . after `Vortex.Point`, which is needed to create
an array of point vortices
=#
v = Vortex.Point.(zvort,Γvort)

# Visualize with the usual plot
ψ .= PotentialFlow.streamfunction(zg,v)
p = plot(ψ,g,color=:black,levels=range(-1,1,length=15),xlim=(-2,2),ylim=(-2,2),xlabel=L"x",ylabel=L"y",title="Streamlines of a pair of vortices",show=true)
plot!(p,v) ## This adds markers for the vortices
add_arrows!(p,v)

#=
Let's evaluate the velocity at the origin, which is a convenient point on the
vertical line in the middle. (Again, the last argument is time, which is irrelevant.)
=#
induce_velocity(0.0+im*0.0,v,0.0)

#=
The real part ($u$) is zero (as expected on this line of symmetry) and the imaginary part
 ($v$) is positive, indicating an **upward flow** between the two vortices.
=#

#=
### Another combination: several vortices
"The tools provided in `PotentialFlow` do not require us to explicitly add the flows;
it does it for us. This is especially useful with a lot of vortices. Let's try to make
a circular **vortex patch**. This is formed from a collection of point vortices arranged
in concentric rings. Each vortex will have the same strength.
=#

# Create one patch at $(-1,0)$ with strength 1
xcent = -1.0
ycent = 0.0
Γ = 1.0
R = 0.5
nring = 5

v1 = vortex_patch(xcent,ycent,Γ,R,nring);

# and another patch at $(1,0)$ with strength -1
xcent = 1.0
ycent = 0.0
Γ = -1.0
R = 0.5
nring = 5

v2 = vortex_patch(xcent,ycent,Γ,R,nring);

# We can combine them easily into something called a **tuple**
vortex_system = (v1,v2);

# Now plot them
ψ .= PotentialFlow.streamfunction(zg,vortex_system)
p = plot(ψ,g,color=:black,xlim=(-2,2),ylim=(-2,2),xlabel=L"x",ylabel=L"y",title="Streamlines of two vortex patches",show=true)
plot!(p,vortex_system)
add_arrows!(p,vortex_system)

# We can evaluate the velocity of this system at any point, e.g.,
z_eval = 0.0+0.0*im
induce_velocity(z_eval,vortex_system,0.0)
