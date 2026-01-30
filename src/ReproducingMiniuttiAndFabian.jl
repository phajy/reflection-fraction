using Gradus
using Plots
include("StandardFuncs.jl")

a = 0.998

m = KerrMetric(1.0, a)

d = ThinDisc(0.0, Inf)
h_out = 100
N_h = 30
N = 1000

heights, geods = calc_geods(m, d; N, h_out, N_h)

disc, bh, missed, total, aisco, bisco = count_fractions(geods, Gradus.isco(m))

plot(heights, disc./total, label="Disc")
plot!(heights, aisco./total, label="Above ISCO")
plot!(heights, bisco./total, label="Below ISCO")
plot!(heights, missed./total, label="inf")

xaxis!("Source height (\$r_g\$)", :log10, xlims=(Gradus.isco(m), 100))
yaxis!("Fraction", ylims=(0,1))

plot(heights, aisco./missed, label="Above ISCO")
xaxis!("Source height (\$r_g\$)", :log10, xlims=(1, 100), )
yaxis!("R", ylims=(0,7))
title!("Miniutti and Fabian Fig 5")