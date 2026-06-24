struct funnel_disk{T} <: AbstractThickAccretionDisc{T}
    inner_radius::T
    outer_radius::T
    α::T

end

function Gradus.cross_section(d::funnel_disk, ρ)
    if ρ<d.inner_radius
        return 0
    elseif ρ > d.outer_radius
        return 0
    else
        (ρ - d.inner_radius) * d.α
    end
end

function funnel_disk(
    alpha:: T,
    m::AbstractMetric{T};
    r_inner = nothing,
    r_outer = 100.0
    ) where {T}
    if r_inner === nothing
        r_inner = Gradus.isco(m)
    end
    
    funnel_disk(r_inner, r_outer, alpha)
end


export funnel_disk

