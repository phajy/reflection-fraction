using Gradus

function point_source_geodesics(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model::AbstractCoronaModel;
    callback = Gradus.domain_upper_hemisphere(),
    kwargs...,
)
    """
    Taken from Fergus' thesis code
    """
    solver_kwargs, setup =
        Gradus.EmissivityProfileSetup(Float64, Gradus.PowerLawSpectrum(2.0); kwargs...)
    δs = deg2rad.(range(setup.δmin, setup.δmax, setup.n_samples))
    # we assume a point source
    x, v = Gradus.sample_position_velocity(m, model)
    velfunc = Gradus.polar_angle_to_velfunc(m, x, v, δs)
    gps = tracegeodesics(
        m,
        x,
        velfunc,
        d,
        setup.λmax;
        save_on = false,
        ensemble = Gradus.EnsembleEndpointThreads(),
        callback = callback,
        trajectories = length(δs),
        solver_kwargs...,
    )

    (; angles = δs, gps = gps)
end

function calc_geods(m, d; N, h_out, N_h)
    """
    Taken from Fergus' thesis code
    """
    heights = collect(logrange(Gradus.isco(m) + 1e-2, h_out, N_h))
    geods = map(heights) do h
        @info h
        lp = LampPostModel(; θ = 1e-5, h = h)
        @time point_source_geodesics(m, d, lp; n_samples = N)
    end
    (; heights, geods)
end


function count_fractions(geods, risco)
    """
    Taken from Fergus' thesis code
    """
    hit_counts = map(geods) do geod
        hit_disc = 0
        hit_above_isco = 0
        hit_below_isco = 0
        hit_bh = 0
        missed = 0
        # need to weight them by their angle
        for (ang, g) in zip(geod.angles, geod.gps)

            w = sin(ang)
            if g.status == StatusCodes.IntersectedWithGeometry
                hit_disc += w
                if g.x[2] >= risco
                    hit_above_isco += w
                else
                    hit_below_isco += w
                end
            elseif g.status == StatusCodes.WithinInnerBoundary
                hit_bh += w
            else
                missed += w
            end
        end
        (;
            disc = hit_disc,
            bh = hit_bh,
            missed = missed,
            above_isco = hit_above_isco,
            below_isco = hit_below_isco,
        )
    end

    disc = [i.disc for i in hit_counts]
    bh = [i.bh for i in hit_counts]
    missed = [i.missed for i in hit_counts]
    above_isco = [i.above_isco for i in hit_counts]
    below_isco = [i.below_isco for i in hit_counts]
    (; disc, bh, missed, total = @.(disc + bh + missed), above_isco, below_isco)
end