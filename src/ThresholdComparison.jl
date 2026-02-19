# -------------------------------- #
# Plot Modified Cross Sections
# -------------------------------- #
h_out = 100
N_h = 30
N = 1000

m = KerrMetric(1.0, 0.0)
threshold_vals = [0.01, 0.05, 0.1, 0.5, 1]
plot()

d_st = ShakuraSunyaev(m; eddington_ratio=0.3)
cs_st = Gradus.cross_section.(d_st, r)
plot!(r, cs_st, label="Default", color=:black, ls=:solid)
for th in threshold_vals
    d = ShakuraSunyaev_custom(m; eddington_ratio=0.3, threshold=th)   
    cs = Gradus.cross_section.(d, r)
    plot!(r, cs, label="Threshold = $th", ls=:dash)
end

xaxis!("Radius (r_g)", :log10)
yaxis!("Cross section (r_g)")
title!("Example of modified cross section")
# savefig("data/results/custom_shakura_sunyaev/ModifiedCrossSection.pdf")

# -------------------------------- #  
# Calculate Geods #
# -------------------------------- #
h_out = 500
N_h = 60
N = 1000
m = KerrMetric(1.0, 0.0)

# Thin Disc
d_thin = ThinDisc(0.0, Inf)
heights_thin, geods_thin = calc_geods(m, d_thin; N,h_out, N_h)

# Default Shakura Sunyaev
d_st = ShakuraSunyaev(m; eddington_ratio=0.3)
heights_st, geods_st = calc_geods(m, d_st; N,h_out, N_h)

threshold_vals = [0.01, 0.05, 0.1, 0.5, 1]
threshold_vals = [0.1]
mapped_geods = map(threshold_vals) do th
    d = ShakuraSunyaev_custom(m; eddington_ratio=0.3, threshold=th)
    heights_custom, geods_custom = calc_geods(m, d; N,h_out, N_h)
    (heights_custom, geods_custom)
end

# -------------------------------- #   
# Plot different fractions #
# -------------------------------- #

# Below ISCO
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, bisco./total, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, bisco./total, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], bisco./total, label="Threshold = $(threshold_vals[n]), below ISCO", linestyle=:dash)
end

title!("Below ISCO Fraction")
xaxis!("Radius (r_g)", :log10)
yaxis!("Fraction")
# savefig("data/results/custom_shakura_sunyaev/BISCO_a=$(m.a).pdf")


# Missed / inf
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, missed./total, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, missed./total, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    print(typeof(mg[2]))
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], missed./total, label="Threshold = $(threshold_vals[n])", linestyle=:dash)
end
title!("Missed Fraction")
xaxis!("Radius (r_g)", :log10)
yaxis!("Fraction")
# savefig("data/results/custom_shakura_sunyaev/missed_a=$(m.a).pdf")

# Above ISCO
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

# savefig("data/results/custom_shakura_sunyaev/AISCO_a=$(m.a).pdf")

# R
plot()

disc, bh, missed, total, aisco, bisco = count_fractions(geods_st, Gradus.isco(m))
plot!(heights_st, aisco./missed, label="Default S&S", linestyle=:solid, color=:black)

disc, bh, missed, total, aisco, bisco = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, aisco./missed, label="Thin Disc", linestyle=:solid, color=:saddlebrown)

for (n, mg) in enumerate(mapped_geods)
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], aisco./missed, label="Threshold = $(threshold_vals[n]), below ISCO", linestyle=:dash)
end
title!("R")
xaxis!("Radius (r_g)", :log10)
yaxis!("R")

# savefig("data/results/custom_shakura_sunyaev/R_a=$(m.a).pdf")

# "Correct" R
plot()
disc, bh, missed, total_st, aisco_st, bisco_th = count_fractions(geods_st, Gradus.isco(m))

disc, bh, missed_thin, total_thin, aisco_thin, bisco_thin = count_fractions(geods_thin, Gradus.isco(m))
plot!(heights_thin, aisco_thin./missed_thin, label="Thin Disc", linestyle=:solid)

for (n, mg) in enumerate(mapped_geods)
    disc, bh, missed, total, aisco, bisco = count_fractions(mg[2], Gradus.isco(m))
    plot!(mg[1], aisco_st./missed, label="Threshold = $(threshold_vals[n])", linestyle=:dash)
end

title!("R")
xaxis!("Radius (r_g)", :log10)
yaxis!("R", ylims=(0,7))