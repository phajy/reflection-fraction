using Gradus

struct ShakuraSunyaev_custom{T} <: AbstractThickAccretionDisc{T}
    Ṁ_Ṁedd::T
    inv_η::T
    inner_radius::T
    r_isco::T
    threshold::T
end

function Gradus.cross_section(d::ShakuraSunyaev_custom, ρ)
    if ρ < d.r_isco
        return d.threshold
    end
    3 * d.inv_η * d.Ṁ_Ṁedd * (1 - sqrt(d.r_isco / ρ))
end

function ShakuraSunyaev_custom(
    m::AbstractMetric{T};
    eddington_ratio = 0.3,
    η = nothing,
    contra_rotating = false,
    threshold=0.1
) where {T}
    r_isco = Gradus.isco(m)
    radiative_efficiency = if isnothing(η)
        1 - CircularOrbits.energy(
            m,
            SVector{2}(r_isco, π / 2);
            contra_rotating = contra_rotating,
        )
    else
        η
    end
    ShakuraSunyaev_custom(T(eddington_ratio), inv(radiative_efficiency), 0.0, r_isco, threshold)
end

export ShakuraSunyaev_custom