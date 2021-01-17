## Potential flow routines

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

    return Vortex.Point.(zvort,strength/cnt)

end
