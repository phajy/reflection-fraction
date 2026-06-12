using Gradus, Plots
include("StandardFuncs.jl")
using ColorSchemes
using Interpolations
# -------- #
# Settings #
# -------- #
default(palette = palette(:Dark2_8))
default(palette = palette(:Set3_12))

save_dir = "$(pwd())/data/results/pub/figures" # Set where to save figures
load_dir = "$(pwd())/data/results/pub/stash" # Set where to load data from
# file_pref = "run_2" # Add a prefix to the filename being saved/loaded
file_pref_load = "run_2"
# overwrite = false

# Initialize grid

# Spin values for which results are calculated / loaded from
a_vals = [0.0, 0.3, 0.5, 0.7, 0.900, 0.990, 0.998]

# Eddington ratios for which results are calculated / loaded from
m_edd_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 

# Coronal heights
h_out = 100 # Maximum height of corona
N_h = 40 # Number of coronal heights considered
N = 10000 # Number of photons that are traced

# Max R as a function of m_edd
begin
    plot(xlabel="m_edd", ylabel="R", title="Maximum R as a function of m_edd", legend=:bottomright, legendfontsize=7)
    for a_sel in a_vals
        m = KerrMetric(1.0, a_sel)
        heights = create_heights(m, h_out, N_h)
        heights_fine = create_heights(m, h_out, 500)
        m_edd_full = copy(m_edd_vals)
        pushfirst!(m_edd_full, 0.0)
        
        R_max = zeros(size(m_edd_full))
        for (n,medd_sel) in enumerate(m_edd_full)
            # Load saved fractions then calculate R
            
            # m_R, m_R_h = calc_max_R(d, heights, heights_fine; ISCO=Gradus.isco(m))
            # R_max[n+1] = m_R
            if medd_sel == 0.0
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-0.0-thin")
                print("Here")
            else
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-$medd_sel-SS")
            end
            m_R, m_R_h = calc_max_R(d, heights, heights_fine; ISCO=Gradus.isco(m))
            R_max[n] = m_R                        
        end
        plot!(m_edd_full, R_max, label="a = $a_sel", ls=:dash, markershape=:circle, ms=2.5)    
    end
    yaxis!(ylims=(1,50), :log10, yticks=([1, 2, 5, 10, 20,  50], [1, 2, 5, 10, 20,  50]))
    # savefig("$save_dir/R_vs_medd.pdf")
end

# ------------------------------ #
# Thin disk reflection fractions #
# ------------------------------ #

begin
pl = plot(grid=true, minorgrid=true)
# colors = [:red, :blue, :green, :cyan, :purple, :orange]
for (n, a) in enumerate(a_vals)
    
    m = KerrMetric(1.0, a)
    heights = create_heights(m, h_out, N_h)

    fname = "$file_pref_load-ref-frac-$a-0.0-thin"
    d_thin = loaddata(load_dir, fname)

    plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="a = $a", ls=:solid) 
end
xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20], [1, 2, 5, 10, 20]))
title!("Reflection Fractions")
display(pl)
end

begin
pl = plot(grid=true, minorgrid=true)

colors = [:blue, :black, :maroon, :green, :pink, :purple, :teal, :orange, :red, :silver]
# Set which accretion rate is being plotted
a_sel = 0.998

# plot!([], [], label="Thin Disk", ls=:dash, c=:black)
# plot!([], [], label="M_edd = $m_edd_sel", ls=:solid, c=:black)


for (n, m_edd) in enumerate(m_edd_vals)
    
    m = KerrMetric(1.0, a_sel)
    heights = create_heights(m, h_out, N_h)

    d_custom = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-$m_edd-SS")

    plot!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="m_edd = $m_edd", ls=:solid, c=colors[n])     

    if n == 1
        d_thin = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-0.0-thin")
        plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="m_edd = 0.0", ls=:solid, c=colors[end]) 
    end
    # break
end

xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20, 25, 50], [1, 2, 5, 10, 20, 25, 50]))
title!("Reflection Fractions, \$ a = $a_sel\$")
# savefig("data/results/custom_shakura_sunyaev/Preliminary/R_comparison_medd$m_edd_sel.pdf")
display(pl)
end

# Maximum reflection fractions

begin
    colors = [:red, :blue, :green, :cyan, :purple, :orange]
    # m_edd_sel = 0.6
    Ignore_below_ISCO = false

    max_R = zeros(size(a_vals)) # Stores max R for given spin
    max_R_h = zeros(size(a_vals)) # Stores height at which max R is reached
    ISCO_vals = zeros(size(a_vals))
    EH_vals = zeros(size(a_vals))

    plot()
    for m_edd_sel in m_edd_vals
        for (n, a) in enumerate(a_vals)    
            
            # Define m to create heights
            m = KerrMetric(1.0, a)
            heights = create_heights(m, h_out, N_h)
            heights_fine = create_heights(m, h_out, 500)
            
            ISCO = Gradus.isco(m) # !! Should make ISCO grid more fine, not just the 5 values of a
            ISCO_vals[n] = ISCO
            EH_vals[n] = Gradus.inner_radius(m) +1e-2

            # Load saved fractions then calculate R
            d_custom = loaddata(load_dir, "$file_pref_load-ref-frac-$a-$m_edd_sel-SS")
            R = d_custom["above_isco"] ./ d_custom["missed"]

            # Convert all NaN values of R into 0s
            # NaNs occur when no photons reach infinity (i.e. very low heights)
            R[isnan.(R)] .= 0
            interp = interpolate((heights,), R, Gridded(Linear()))
            R_interpolated = interp.(heights_fine)
            # Only consider heights greater than ISCO
            
            h_above_ISCO = heights_fine[heights_fine .>= ISCO]
            R_above_ISCO = R_interpolated[heights_fine .> ISCO]

            # Calcualte max values
            if !Ignore_below_ISCO
                max_R[n] = findmax(R_interpolated)[1]
                max_R_h[n] = heights_fine[findmax(R_interpolated)[2]]
            else
                max_R[n] = findmax(R_above_ISCO)[1]
                max_R_h[n] = h_above_ISCO[findmax(R_above_ISCO)[2]]
            end
        end
        scatter!(a_vals, max_R_h, label="m_edd = $m_edd_sel", shape=:rect, ms=2, alpha=.7)
    end
    
    # scatter!(a_vals, max_R, label="Thick disk", shape=:rect, ms=5, alpha=.7)

    for (n, a) in enumerate(a_vals)    
        # Define m to create heights
        m = KerrMetric(1.0, a)
        heights = create_heights(m, h_out, N_h)
        ISCO = Gradus.isco(m)

        # Load saved fractions then calculate R
        d_thin = loaddata(load_dir, "$file_pref_load-ref-frac-$a-0.0-thin")
        R = d_thin["above_isco"] ./ d_thin["missed"]

        # Convert all NaN values of R into 0s
        # NaNs occur when no photons reach infinity (i.e. very low heights)
        R[isnan.(R)] .= 0

        # Only consider heights greater than ISCO
        h_above_ISCO = heights[heights .>= ISCO]
        R_above_ISCO = R[heights .> ISCO]
        
        # scatter!(h_above_ISCO, R_above_ISCO, shape=:rect)
        # plot!(heights, R)

        # Calcualte max values
        if !Ignore_below_ISCO
            max_R[n] = findmax(R)[1]
            max_R_h[n] = heights[findmax(R)[2]]
        else
            max_R[n] = findmax(R_above_ISCO)[1]
            max_R_h[n] = h_above_ISCO[findmax(R_above_ISCO)[2]]
        end
    end

    # scatter!(a_vals, max_R, label="Thin Disk")
    
    scatter!(a_vals, max_R_h, label="Thin Disk")
    plot!(a_vals, ISCO_vals, label="ISCO", ls=:dash, c=:black)
    plot!(a_vals, EH_vals, label="EH", ls=:dash, c=:black)

    xaxis!("Spin")
    yaxis!("Source Height (r_g)")
    title!("Source height at Maximum R")
    # savefig("data/results/custom_shakura_sunyaev/Preliminary/R_max_h_medd$m_edd_sel-IgnoreISCO-$Ignore_below_ISCO.pdf")
end

yaxis!(:log10, ylims=(1, 4))

# Interpolation
using Interpolations


function calc_max_R(data, heights, heights_fine; ISCO = -1.0)
    R = data["above_isco"] ./ data["missed"]

    # Convert all NaN values of R into 0s
    # NaNs occur when no photons reach infinity (i.e. very low heights)
    R[isnan.(R)] .= 0
    interp = interpolate((heights,), R, Gridded(Linear()))
    R_interpolated = interp.(heights_fine)

    # Calcualte max values
    if ISCO == -1
        max_R = findmax(R_interpolated)[1]
        max_R_h = heights_fine[findmax(R_interpolated)[2]]
    else
        h_above_ISCO = heights_fine[heights_fine .>= ISCO]
        R_above_ISCO = R_interpolated[heights_fine .> ISCO]
        max_R = findmax(R_above_ISCO)[1]
        max_R_h = h_above_ISCO[findmax(R_above_ISCO)[2]]
    end
    return max_R, max_R_h
end






begin
pl = plot(grid=true, minorgrid=true)

colors = [:blue, :black, :maroon, :green, :pink, :purple, :teal, :orange, :red, :silver]
# Set which accretion rate is being plotted
m_edd_sel = 1.0

plot!([], [], label="Thin Disk", ls=:dash, c=:black)
plot!([], [], label="M_edd = $m_edd_sel", ls=:solid, c=:black)

for (n, a) in enumerate(a_vals)
    
    m = KerrMetric(1.0, a)
    heights = create_heights(m, h_out, N_h)

    d_custom = loaddata(load_dir, "$file_pref_load-ref-frac-$a-$m_edd_sel-SS")

    plot!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="a = $a", ls=:solid, c=colors[n]) 

    d_thin = loaddata(load_dir, "$file_pref_load-ref-frac-$a-0.0-thin")

    plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="", ls=:dash, c=colors[n]) 
    # break
end

xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20, 25], [1, 2, 5, 10, 20, 25]))
title!("Reflection Fractions, \$\\dot{m} = $m_edd_sel\$")
# savefig("data/results/custom_shakura_sunyaev/Preliminary/R_comparison_medd$m_edd_sel.pdf")
display(pl)
end