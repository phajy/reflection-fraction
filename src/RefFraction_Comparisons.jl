using Gradus
using Plots

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


#### -------------------------------------------------- ####
#### Comparing Thin Disk to Thick Disks, maximal spin   ####
#### -------------------------------------------------- ####

m = KerrMetric(1.0, 0.998)
d = ThinDisc()
h_out = 100
N_h = 30
N = 1000

heights_1, geods_1 = calc_geods(m, d; N, h_out, N_h)

m = KerrMetric(1.0, 0.998)
d = ShakuraSunyaev(m; eddington_ratio=.1)
h_out = 100
N_h = 30
N = 1000
heights_2, geods_2 = calc_geods(m, d; N,h_out, N_h)

m = KerrMetric(1.0, 0.998)
d = ShakuraSunyaev(m; eddington_ratio=.3)
h_out = 100
N_h = 30
N = 1000
heights_3, geods_3 = calc_geods(m, d; N,h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_1, Gradus.isco(m))
plot(heights_1, aisco./total, label="Thin Disk, above ISCO", linestyle=:solid, color=:green, legend=:right)
plot!(heights_1, bisco./total, label="Thin Disk, below ISCO", linestyle=:dash, color=:green)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_2, Gradus.isco(m))
plot!(heights_2, aisco./total, label="Mdot = 0.1, above ISCO", linestyle=:solid, color=:red)
plot!(heights_2, bisco./total, label="Mdot = 0.1, below ISCO", linestyle=:dash, color=:red)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_3, Gradus.isco(m))
plot!(heights_2, aisco./total, label="Mdot = 0.3, above ISCO", linestyle=:solid, color=:blue)
plot!(heights_2, bisco./total, label="Mdot = 0.3, below ISCO", linestyle=:dash, color=:blue)

m = KerrMetric(1.0, 0.998)
vline!([Gradus.isco(m)], label="ISCO", color=:black, ls=:dot)

title!("a = 0.998")
xaxis!("Source Height (r_g)", :log10)
ylabel!("Fraction (photon count / total)")

#### -------------------------------------------------- ####
### Comparing Thin Disk to Thick Disks, 0 spin ###
#### -------------------------------------------------- ####

m = KerrMetric(1.0, 0.)

d = ThinDisc()
h_out = 100
N_h = 30
N = 1000
heights_1, geods_1 = calc_geods(m, d; N, h_out, N_h)

d = ShakuraSunyaev(m; eddington_ratio=.1)
h_out = 100
N_h = 30
N = 1000
heights_2, geods_2 = calc_geods(m, d; N,h_out, N_h)

# m = KerrMetric(1.0, 0.998)
d = ShakuraSunyaev(m; eddington_ratio=.3)
h_out = 100
N_h = 30
N = 1000
heights_3, geods_3 = calc_geods(m, d; N,h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_1, Gradus.isco(m))
plot(heights_1, aisco./total, label="Thin Disk, above ISCO", linestyle=:solid, color=:green, legend=:right)
plot!(heights_1, bisco./total, label="Thin Disk, below ISCO", linestyle=:dash, color=:green)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_2, Gradus.isco(m))
plot!(heights_2, aisco./total, label="Mdot = 0.1, above ISCO", linestyle=:solid, color=:red)
plot!(heights_2, bisco./total, label="Mdot = 0.1, below ISCO", linestyle=:dash, color=:red)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_3, Gradus.isco(m))
plot!(heights_2, aisco./total, label="Mdot = 0.3, above ISCO", linestyle=:solid, color=:blue)
plot!(heights_2, bisco./total, label="Mdot = 0.3, below ISCO", linestyle=:dash, color=:blue)

# m = KerrMetric(1.0, 0.998)
vline!([Gradus.isco(m)], label="ISCO", color=:black, ls=:dot)

# m = KerrMetric(1.0, 0.0)
# vline!([Gradus.isco(m)], label="ISCO, a = 0.0")

title!("a = 0.")
xaxis!("Source Height (r_g)", :log10)
ylabel!("Fraction (photon count / total)")

#### -------------------------------------------------- ####
### Comparing large vs small spins ###
#### -------------------------------------------------- ####

m = KerrMetric(1.0, 0.998)
d = ShakuraSunyaev(m; eddington_ratio=.8)
h_out = 100
N_h = 30
N = 1000

heights_1, geods_1 = calc_geods(m, d; N, h_out, N_h)

m = KerrMetric(1.0, 0.)
d = ShakuraSunyaev(m; eddington_ratio=.8)
h_out = 100
N_h = 30
N = 1000

heights_2, geods_2 = calc_geods(m, d; N,h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_1, Gradus.isco(m))
plot(heights_1, aisco./total, label="a=0.998, above ISCO", linestyle=:dash)
plot!(heights_1, bisco./total, label="a=0.998, below ISCO", linestyle=:dash)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_2, Gradus.isco(m))
plot!(heights_2, aisco./total, label="a=0., above ISCO", linestyle=:dash)
plot!(heights_2, bisco./total, label="a=0., below ISCO", linestyle=:dash)

m = KerrMetric(1.0, 0.998)
vline!([Gradus.isco(m)], label="ISCO, a = 0.998")

m = KerrMetric(1.0, 0.0)
vline!([Gradus.isco(m)], label="ISCO, a = 0.0")

title!("Eddington Ratio = 0.3")
xaxis!("Heights (r_g)", :log10)
ylabel!("Fraction (photon count / total)")

#### -------------------------------------------------- ####
#### Going from thin disc to thick disk - WIP!!
#### -------------------------------------------------- ####

Mdot_vals = [.1, 0.2, 0.3, 0.5, 0.7, 0.8, .9]

m = KerrMetric(1.0, 0.0)

h_out = 100
N_h = 30
N = 1000

a = plot(legend=:topright)
for mdot in Mdot_vals
    d = ShakuraSunyaev(m; eddington_ratio=mdot)
    heights, geods = calc_geods(m, d; N,h_out, N_h)

    disc, bh, missed, total, aisco, bisco = count_fractions(geods, Gradus.isco(m))
    plot!(a, heights, aisco./total, label="Mdot = $mdot", linestyle=:dash)
    # plot!(heights_1, bisco./total, label="a=0.998, below ISCO", linestyle=:dash)
end

title!(a, "a = 0.0")
xaxis!("Source height (\$r_g\$)", :log10, xlims=(Gradus.isco(m),100))
yaxis!("Fraction", ylims=(0,1))
vline!(a, [Gradus.isco(m)], label="ISCO, a = 0.0", color=:black)
hline!([0.5], ls=:dot, c=:black)
display(a)