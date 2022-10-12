# Define structures and functions for extending ViscousFlow functions
# to work on PotentialFlow elements (just for plotting and basic calculation)

export GridCache

import ViscousFlow: Nodes, Edges, Primal, Dual, BasicILMCache, x_gridcurl, y_gridcurl,
                    x_gridgrad, y_gridgrad, PhysicalGrid, SurfaceScalarCache,
                    zeros_gridcurl, zeros_gridgrad

PotFlowElements = Union{PotentialFlow.Element,Vector{V} where {V<:PotentialFlow.Element},Tuple}

GridCache(g::PhysicalGrid) = SurfaceScalarCache(g)

function ViscousFlow.streamfunction(v::PotFlowElements,cache::BasicILMCache;angle=nothing)
    zg = x_gridcurl(cache) + im*y_gridcurl(cache)
    if !isnothing(angle)
        rot = exp(-im*(π+angle))  # rotation operator, which moves rotates from $\tau$ to $-\pi$.
        zg .*= rot
    end
    ψ = zeros_gridcurl(cache)
    ψ .= PotentialFlow.streamfunction(zg,v)
end

function ViscousFlow.streamfunction!(ψ::Nodes{Dual},v::PotFlowElements,cache::BasicILMCache;kwargs...)
  ψ .= ViscousFlow.streamfunction(v,cache;kwargs...)
end

function ViscousFlow.velocity(v::PotFlowElements,cache::BasicILMCache)
    x = x_gridgrad(cache)
    y = y_gridgrad(cache)
    zgu = x.u + im*y.u
    zgv = x.v + im*y.v
    vel = zeros_gridgrad(cache)
    vel.u .= real(induce_velocity(zgu,v,0.0))
    vel.v .= imag(induce_velocity(zgv,v,0.0))
    return vel
end

function ViscousFlow.velocity!(vel::Edges{Primal},v::PotFlowElements,cache::BasicILMCache;kwargs...)
  vel .= ViscousFlow.velocity(v,cache;kwargs...)
end
