using Gradus
using Plots
include("StandardFuncs.jl")
include("CustomShakuraSunyaev.jl")

# -------------------------------- #
# Calculate results for 3 values of spin
# -------------------------------- #

a_vals = [0.0, 0.5, 0.998]

h_out = 100
N_h = 30
N = 1000
thrshld = 0.1
m_edd = 0.3

# Calculate and store values to file

for a in a_vals
    m = KerrMetric(1.0, a)
    d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
    heights, geods = calc_geods(m, d; N,h_out, N_h)
    cf = count_fractions(geods, Gradus.isco(m))
    stashdata("ref-frac-$a-$m_edd-custom"; cf...)

    d_thin = ThinDisc(0.0, Inf)
    heights_thin, geods_thin = calc_geods(m, d_thin; N,h_out, N_h)
    cf = count_fractions(geods_thin, Gradus.isco(m))
    stashdata("ref-frac-$a-0.0-thin"; cf...)

    # Default Shakura Sunyaev
    d_st = ShakuraSunyaev(m; eddington_ratio=m_edd)
    heights_st, geods_st = calc_geods(m, d_st; N,h_out, N_h)
    cf = count_fractions(geods_st, Gradus.isco(m))
    stashdata("ref-frac-$a-$m_edd-default"; cf...)
end
plot()
# Read from files and plot
begin
c = plot(legend=:right)
for str in ["disc", "bh", "missed"]
    # str = "disc"
    for a in a_vals[1:1]
        m = KerrMetric(1.0, a)
        heights = collect(logrange(Gradus.isco(m) + 1e-2, h_out, N_h))

        d_custom = loaddata("ref-frac-$a-$m_edd-custom")
        d_thin = loaddata("ref-frac-$a-0.0-thin")
        d_default = loaddata("ref-frac-$a-$m_edd-default")
        
        R_thin = d_thin["above_isco"] ./ d_thin["missed"]
        R_SS = (d_default["above_isco"] ./ d_default["total"]) ./ (d_custom["missed"] ./ d_custom["total"])

        # plot!(heights, R_thin, label="Thin Disc")
        # plot!(heights, R_SS, label="S&S")

        plot!(heights, d_custom[str] ./ d_custom["total"], label="$str, custom", ls=:dash)
        plot!(heights, d_thin[str] ./ d_thin["total"] ,label="$str, thin")
        # plot!(heights, d_default[str] ./ d_default["total"] ,label="default")

        # plot!(heights, d_custom["total"], label="custom")
        # plot!(heights, d_thin["total"], label="thin")
        # plot!(heights, d_default["total"], label="default")

        # plot!(heights, d_default["above_isco"] ./ d_default["total"], label="Above a=$a")
        # plot!(heights, d_thin["above_isco"] ./ d_thin["total"] ,label="Thin above a=$a")
        
        title!("Fractions, a = $a")
        
        # break
        
    end
xaxis!("Source Height (r_g)", :log10, minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction")
display(c)
# savefig("FractionsBug_ExampleFig.png")
end


end

### Comparing with Fergus' definition of disc ###

function FilledDisc(m)
    ss = ShakuraSunyaev(m)
    r_isco = Gradus.isco(m)
    ThickDisc() do x
        if x > r_isco
            Gradus.cross_section(ss, x)
        else
            0.01
        end
    end
end

a = 0.0
h_out = 100
N_h = 30
N = 5000
thrshld = 0.01
m_edd = 0.1

m = KerrMetric(1.0, a)

d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
heights, geods_custom = calc_geods(m, d; N,h_out, N_h)
frac_custom = count_fractions(geods_custom, Gradus.isco(m))

d = ThinDisc(0.0, 500.0)
heights, geods_thin = calc_geods(m, d; N,h_out, N_h)
frac_thin = count_fractions(geods_thin, Gradus.isco(m))

# d = FilledDisc(m)
# heights, geods = calc_geods(m, d; N,h_out, N_h)
# frac_filled = count_fractions(geods, Gradus.isco(m))


plot(legend=:right)
plot!(heights, frac_custom.missed ./ frac_custom.total, ls=:dash, label="Infinity")
plot!(heights, frac_thin.missed ./ frac_thin.total, ls=:solid, label="")
# plot!(heights, frac_filled.missed ./ frac_filled.total, ls=:dashdot, label="")

plot!(heights, frac_custom.disc ./ frac_custom.total, ls=:dash, label="Disc",)
plot!(heights, frac_thin.disc ./ frac_thin.total, ls=:solid, label="",)
# plot!(heights, frac_filled.disc ./ frac_filled.total, ls=:dashdot, label="",)

plot!(heights, frac_custom.bh ./ frac_custom.total, ls=:dash, label="Above ISCO",)
plot!(heights, frac_thin.bh./ frac_thin.total, ls=:solid, label="",)
# plot!(heights, frac_filled.bh ./ frac_filled.total, ls=:dashdot, label="",)

xaxis!("Source Height (r_g)", :log10, minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction")
title!("Eddington ratio = $m_edd, a=$a")
savefig("Fractions_EddRatio0.1.pdf")

### Plot Geometry ###
plot()
function circleShape(h, k, r)
    θ = LinRange(0, 2*π, 500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

plot(circleShape(0,0,Gradus.event_horizon(m)[1][1]), seriestype=[:shape,], aspect_ratio=1, 
c=:black, label="")
plot!(circleShape(0,10,.1), seriestype=[:shape,], aspect_ratio=1, c=:blue, label="")

r = LinRange(0, 100, 500)
d_st = ShakuraSunyaev(m; eddington_ratio=0.3)
cs_st = Gradus.cross_section.(d_st, r)
plot!(r, cs_st, label="Default", color=:black, ls=:solid)

d = ShakuraSunyaev_custom(m; eddington_ratio=0.3, threshold=.1)   
cs = Gradus.cross_section.(d, r)
plot!(r, cs, label="Threshold = $thrshld", ls=:dash)

x_vec_custom = zeros(size(geods_custom[30].gps))
for (i, gp) in enumerate(geods_custom[1].gps)
    x_vec_custom[i] = gp.x[2]
    # println(typeof(gp.x[2]))
end
x_vec_thin = zeros(size(geods_thin[30].gps))
for (i, gp) in enumerate(geods_thin[1].gps)
    x_vec_thin[i] = gp.x[2]
    # println(typeof(gp.x[2]))
end
x_vec
histogram(x_vec_thin, bins=logrange(0.01, 10000, 100), alpha=1, label="thin")
histogram!(x_vec_custom, bins=logrange(0.01, 10000, 100), alpha=0.8, label="Custom")
xaxis!(xlims=(8e3, 10000), :log10)
# yaxis!(ylims=(.1, .4))
# (geods[3].gps[1].x[2])
typeof(geods[3].gps)


### ----------------------------------------- ###
##### Plot paths #####
### ----------------------------------------- ###

#This requires the tracegeodesics function to be modified such that it saves the traced path
#Also requires removing the ensemblethreads... solving option.

begin
a = 0.0
h_out = 100
N_h = 2
N = 100
thrshld = 0.01
m_edd = 0.3

m = KerrMetric(1.0, a)

d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
heights, geods_custom = calc_geods(m, d; N,h_out, N_h)

d = ThinDisc(0.0, 500.0)
heights, geods_thin = calc_geods(m, d; N,h_out, N_h)

print(":)") # To ensure that it doesn't print out all the returned geodesics :)
end

# The counting fucntion doesn't directly work when we return the entire traced geodesic #
# frac_custom = count_fractions(geods_custom, Gradus.isco(m))
# frac_thin = count_fractions(geods_thin, Gradus.isco(m))


pl_custom = plot_paths(geods_custom[2].gps, legend=false, extent=200, N_points=5000, plane=:XZ, t_span=4000)
title!("Traced Paths, Custom S&S disk")
xaxis!(xlims=(-600, 0), xlabel="Radius (rg)")
yaxis!(ylims=(-10, 200), ylabel="Height (rg)")
hline!([0], c=:black, ls=:dash)
vline!([-500], ls=:dash, c="black")

pl_thin = plot_paths(geods_thin[2].gps, legend=false, extent=200, N_points=6000, plane=:XZ, t_span=4000)

title!("Traced Paths, Thin Disk")
xaxis!(xlims=(-600, 0), xlabel="Radius (rg)")
yaxis!(ylims=(-10, 200), ylabel="Height (rg)")
hline!([0], c=:black, ls=:dash)