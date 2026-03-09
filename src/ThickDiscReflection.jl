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
N = 10_000
thrshld = 0.1
m_edd = 0.6

# Calculate and store values to file
for a in a_vals
    m = KerrMetric(1.0, a)
    # d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
    d = ShakuraSunyaev(m; eddington_ratio=m_edd) ∘ ThinDisc(0.0, Inf)
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


# Plot reflection fractions #
begin
plot()
colors = [:red, :blue, :green]
a = a_vals[3]
for (n, m_edd) in enumerate([0.3, 0.6])
    
    m = KerrMetric(1.0, a)
    heights = collect(logrange(Gradus.isco(m) + 1e-2, h_out, N_h))
    
    d_custom = loaddata("ref-frac-$a-$m_edd-custom")
    d_thin = loaddata("ref-frac-$a-0.0-thin")
    d_default = loaddata("ref-frac-$a-$m_edd-default")

    if n == 1
        plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], ls=:dash, c=:black, label="Thin Disk")
    end

    scatter!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="M_edd = $m_edd", ls=:solid, c=colors[n]) 
    # plot!(heights, d_default["above_isco"] ./ d_custom["missed"], label="R, custom 2") 
    
end

xaxis!("Source Height (r_g)", :log10,minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction",)
title!("Reflection Fractions, a = $a")
savefig("data/results/custom_shakura_sunyaev/Preliminary/R_CompositeModel_a$a.pdf")
end


# Plot photon fractions #
plot()
begin
c = plot(legend=:topright)
for str in ["missed", "above_isco", "disc"]
    # str = "disc"
    for a in a_vals[1:1]
        m = KerrMetric(1.0, a)
        heights = collect(logrange(Gradus.isco(m) + 1e-2, h_out, N_h))

        d_custom = loaddata("ref-frac-$a-$m_edd-custom")
        d_thin = loaddata("ref-frac-$a-0.0-thin")
        d_default = loaddata("ref-frac-$a-$m_edd-default")
        
        # R_thin = d_thin["above_isco"] ./ d_thin["missed"]
        # R_SS = (d_default["above_isco"] ./ d_default["total"]) ./ (d_custom["missed"] ./ d_custom["total"])

        # plot!(heights, R_thin, label="Thin Disc")
        # plot!(heights, R_SS, label="S&S")

        # plot!(heights, d_custom[str] ./ d_custom["total"], label="$str, custom", ls=:dash)
        # plot!(heights, d_thin[str] ./ d_thin["total"] ,label="$str, thin")
        # plot!(heights, d_default[str] ./ d_default["total"] ,label="default")

        # plot!(heights, d_custom["total"], label="custom")
        # plot!(heights, d_thin["total"], label="thin")
        # plot!(heights, d_default["total"], label="default")

        # plot!(heights, d_default["above_isco"] ./ d_default["total"], label="Above a=$a")
        # plot!(heights, d_thin["above_isco"] ./ d_thin["total"] ,label="Thin above a=$a")
        
        # # R
        plot!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="R, custom") 
        # plot!(heights, d_default["above_isco"] ./ d_custom["missed"], label="R, custom 2") 
        plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="R, thin")
        
        # Ratios
        # plot!(heights, d_custom[str] ./ d_thin[str], label="$str, custom-thin", ls=:dash)
        # plot!(heights, d_default[str] ./ d_custom[str], label="$str, default-custom", ls=:dash)

        title!("Fractions, a = $a, m_edd = $m_edd")
    end
xaxis!("Source Height (r_g)", :log10,minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction",)
# xlims!(1, 100)
# ylims!(0.9, 1.1)
display(c)
break
# savefig("FractionsBug_ExampleFig.png")
end
hline!([1], c=:black, ls=:dot)
end




### Quicker code for testing purposes (samples fewer photons) ###

a = 0.0
h_out = 100
N_h = 30
N = 1000
thrshld = 0.1
m_edd = 0.3

m = KerrMetric(1.0, a)

composite_model = ShakuraSunyaev(m) ∘ ThinDisc(0.0, Inf)

# d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
d = composite_model
heights, geods_custom = calc_geods(m, d; N,h_out, N_h)
frac_custom = count_fractions(geods_custom, Gradus.isco(m))

d = ThinDisc(0.0, Inf)
heights, geods_thin = calc_geods(m, d; N,h_out, N_h)
frac_thin = count_fractions(geods_thin, Gradus.isco(m))

d = ShakuraSunyaev(m)
heights, geods_default = calc_geods(m, d; N,h_out, N_h)
frac_default = count_fractions(geods_thin, Gradus.isco(m))

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

plot!(heights, frac_custom.bh ./ frac_custom.total, ls=:dash, label="BH",)
plot!(heights, frac_thin.bh./ frac_thin.total, ls=:solid, label="",)
# plot!(heights, frac_filled.bh ./ frac_filled.total, ls=:dashdot, label="",)

plot(heights, frac_custom.above_isco ./ frac_custom.missed,)
plot!(heights, frac_thin.above_isco ./ frac_thin.missed,)

xaxis!("Source Height (r_g)", :log10, minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction")
title!("Eddington ratio = $m_edd, a=$a")
# savefig("Fractions_EddRatio0.1.pdf")

# -------------------------------- #
# -------------------------------- #
### Plot Geometry                ###
# -------------------------------- #
# -------------------------------- #

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
