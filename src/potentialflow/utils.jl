## Potential flow routines
using LinearAlgebra: diagm


function complexgrid(x::AbstractVector,y::AbstractVector)
    z = zeros(ComplexF64,length(x),length(y))
    @. z = x + im*y'
    return z
end

"""
    vortex_patch(xcent,ycent,strength,radius,nring) -> Vector{Vortex.Point}

Create a list of point vortices in the form of a vortex patch, a set of `nring` concentric
rings centered at `(xcent,ycent)` with radius `radius`. The overall circulation of the
patch is defined by `strength`.
"""
function vortex_patch(xcent,ycent,strength,radius,nring)
    Δr = radius/(nring-1)
    zcent = xcent+im*ycent

    r = 0.0

    zvort = ComplexF64[]
    cnt = 0
    for i = 1:nring
        θ = 0.0
        nv = max(1,8*(i-1))
        r = (i-1)*Δr
        for j = 1:nv
            push!(zvort,zcent + r*exp(im*θ))
            cnt += 1
            θ += 2π/nv
        end
    end

    return PotentialFlow.Vortex.Point.(zvort,strength/cnt)

end


"""
    dotproduct(a,b)

Compute the dot product between `a` and `b`, each given as complex
forms of two-dimensional vectors.
"""
dotproduct(a::Number, b::Number) = real(conj(a)*b)


# setup and solution of system of equations
function simulate_flow(unit_sources, Δslist, n̂, other_elements; tracer_start = collect(-3.0 .+ range(-3,3,length=31)*im), Tmax = 20.0, Δt = 0.01)

    # locations at which we enforce no flow through
    targets = PotentialFlow.Elements.position.(unit_sources)

    # velocity induced by other elements at the target points
    other_vel = induce_velocity(targets,other_elements,0)

    b = -dotproduct.(n̂, other_vel)

    A = [dotproduct(n, induce_velocity(target, unit_source, 0)) for (n, target) in zip(n̂, targets), unit_source in unit_sources]
    A .+= 0.5*diagm(0 => 1 ./ Δslist)

    Q = A \ b

    sources = PotentialFlow.Source.Point.(targets, Q)

    # compute other stuff, like surface velocity, pressure, and tracer trajectories

    us = surface_velocity(targets, sources, other_elements, n̂, Δslist)

    Cp = 1 .- dotproduct.(us,us) #/dotproduct(U∞,U∞)

    tx, ty = compute_trajectories((other_elements,sources),tracer_start,Tmax,Δt=Δt)


    sources, us, Cp, tx, ty
end

function surface_velocity(targets,sources,other_elements,n̂,Δslist)
    us = induce_velocity(targets, (other_elements,sources), 0)
    for (i,(ni,Δsi,source)) in enumerate(zip(n̂,Δslist,sources))
        us[i] -= 0.5*im*source.S*ni/Δsi
    end
    return us
end

# some shapes

function equilateraltriangle(center::Number,Δs_nom::Real;f=1,len::Real=2)
    half = 0.5*len
    shift = Int(2*f)
    n = ceil(2*half/Δs_nom)-shift
    Δs_actual = float(2*half/(n+shift))
    s = -(half - f*Δs_actual):Δs_actual:(half - f*Δs_actual)

    centershift = center - im*√3/3

    bottom = copy(s) .+ centershift
    right = (s .+ half).*exp(im*2π/3) .+ half .+ centershift
    left = (s .+ half)*exp(-im*2π/3) .+ im*√3*half .+ centershift


    z = vcat(bottom,right,left)

    n̂ = vcat(
        fill(-1.0im, length(bottom)),
        fill(exp(im*π/6), length(right)),
        fill(exp(im*5π/6), length(left)))

    Δs = fill(Δs_actual, length(n̂))

    return z, n̂, Δs
end

function circle(center::Number,radius::Float64,Δs::Real)
    n = ceil(Int,2π*radius/Δs)
    θ = range(0, 2π, length=n+1)[1:end-1]
    z = center .+ radius*exp.(im.*θ)
    n̂ = radius*exp.(im.*θ)
    Δs = radius*fill(step(θ), length(θ))
    return z, n̂, Δs
end
