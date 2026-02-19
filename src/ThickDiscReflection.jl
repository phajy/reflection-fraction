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
begin
plot()
str = "disc"
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

    plot!(heights, d_custom[str] ./ d_custom["total"], label="custom")
    plot!(heights, d_thin[str] ./ d_thin["total"] ,label="thin")
    plot!(heights, d_default[str] ./ d_default["total"] ,label="default")

    # plot!(heights, d_custom["total"], label="custom")
    # plot!(heights, d_thin["total"], label="thin")
    # plot!(heights, d_default["total"], label="default")

    # plot!(heights, d_default["above_isco"] ./ d_default["total"], label="Above a=$a")
    # plot!(heights, d_thin["above_isco"] ./ d_thin["total"] ,label="Thin above a=$a")
    
    title!("$str, a=$a")
    
    # break
    
end

xaxis!("Source Height (r_g)", :log10, minorgrid=true, minroticks=5, xticks=([1, 10, 30, 100], [1, 10, 30, 100]))
yaxis!("Fraction")

end