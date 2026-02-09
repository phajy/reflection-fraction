using Gradus
using Plots
include("StandardFuncs.jl")
include("CustomShakuraSunyaev.jl")

### Hack for extending disc down to the event horizon
h_out = 100
N_h = 60
N = 1000

m = KerrMetric(1.0, 0.998)
plot()
d_st = ShakuraSunyaev(m; eddington_ratio=0.3)
cs_st = Gradus.cross_section.(d_st, r)
plot!(r, cs_st, label="Default", color=:black, ls=:solid)
for th in threshold_vals

    d = ShakuraSunyaev_custom(m; eddington_ratio=0.3, threshold=th)   

    r = LinRange(1, 100, 1000)
    cs = Gradus.cross_section.(d, r)
    plot!(r, cs, label="Threshold = $th", ls=:dash)
end

xaxis!("Radius (r_g)", :log10)
yaxis!("Cross section (r_g)")
title!("Example of modified cross section")
Gradus.isco(m)

### Exploring different threshold values
h_out = 100
N_h = 60
N = 1000
m = KerrMetric(1.0, 0.998)

# Thin Disc
d_thin = ThinDisc(0.0, Inf)
heights_thin, geods_thin = calc_geods(m, d_thin; N,h_out, N_h)

# Default Shakura Sunyaev
d_st = ShakuraSunyaev(m; eddington_ratio=0.3)
heights_st, geods_st = calc_geods(m, d_st; N,h_out, N_h)


# Custom Shakura Sunyaev
threshold_vals = [0.01, 0.05, 0.1, 0.5, 1]

mapped_geods = map(threshold_vals) do th
    d = ShakuraSunyaev_custom(m; eddington_ratio=0.3, threshold=th)
    heights_custom, geods_custom = calc_geods(m, d; N,h_out, N_h)
    (heights_custom, geods_custom)
end

# Plot below ISCO fraction
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, bisco./total, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, bisco./total, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    print(typeof(mg[2]))
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], bisco./total, label="Threshold = $(threshold_vals[n]), below ISCO", linestyle=:dash)
end
title!("Below ISCO Fraction")
xaxis!("Radius (r_g)", :log10)
yaxis!("Fraction")

# Plot missed fraction
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, missed./total, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, missed./total, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    print(typeof(mg[2]))
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], missed./total, label="Threshold = $(threshold_vals[n]), below ISCO", linestyle=:dash)
end
title!("Missed Fraction")
xaxis!("Radius (r_g)", :log10)
yaxis!("Fraction")

# Above ISCO fraction
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, aisco./total, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, aisco./total, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    # print(typeof(mg[2]))
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], aisco./total, label="Threshold = $(threshold_vals[n])", linestyle=:dash)
end
title!("Above ISCO Fraction")
xaxis!("Radius (r_g)", :log10)
yaxis!("Fraction")

# R
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, aisco./missed, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, aisco./missed, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    # print(typeof(mg[2]))
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], aisco./missed, label="Threshold = $(threshold_vals[n]), below ISCO", linestyle=:dash)
end
title!("R")
xaxis!("Radius (r_g)", :log10)
yaxis!("R")



### Fix a value for the threshold = 0.1 and compute R
m = KerrMetric(1.0, 0.0)
plot()

d_thin = ThinDisc(0.0, Inf)
heights_thin, geods_thin = calc_geods(m, d_thin; N,h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, aisco./missed, label="Thin Disc", linestyle=:dash, color=:green)

d = ShakuraSunyaev_custom(m; eddington_ratio=0.1)
heights_custom, geods_custom = calc_geods(m, d; N,h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_custom, Gradus.isco(m))
plot!(heights_custom, aisco./missed, label="Custom S&S, mdot = 0.11", linestyle=:dash, color=:green)

# title!("Comparing below ISCO and missed fractions; a = 0.5")
xaxis!("Radius (r_g)", :log10)
yaxis!("R")
