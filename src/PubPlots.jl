using Gradus, Plots
include("StandardFuncs.jl")
using ColorSchemes
using Interpolations

# todo set it up such that running the entire thing will produce and save all the required figures
# todo color palette

# -------- #
# Settings #
# -------- #

a = default(palette = palette(:Dark2_8))
# colors = [:red, :blue, :green, :cyan, :purple, :orange]
colors = palette(:Dark2_8)
linestyles = [:solid, :dash, :dashdot, :dashdotdot, :dot]

# todo change below to more robust directory trackers
const _ROOT = joinpath(@__DIR__, "..")
save_dir = joinpath(_ROOT, "data", "results", "pub", "figures") # Set where to write data to
load_dir = joinpath(_ROOT, "data", "results", "pub", "stash") # Set where to load data from
# file_pref = "run_2" # Add a prefix to the filename being saved/loaded
file_pref_load = "run_2"
overwrite = true

# Initialize grid #

# Spin values for which results are calculated / loaded from
# a_vals = [0.0, 0.3, 0.5, 0.7, 0.900, 0.990, 0.998]
# a_vals = [0.0, 0.15, 0.3, 0.5, 0.7, 0.8, 0.900, 0.95, 0.990, 0.998]
a_vals = [0.0, 0.1, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.900, 0.95, 0.990, 0.998]

# Eddington ratios for which results are calculated / loaded from
m_edd_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
m_edd_full = copy(m_edd_vals) # Add thin disk case
pushfirst!(m_edd_full, 0.0)

# Coronal heights
h_out = 100 # Maximum height of corona
N_h = 40 # Number of coronal heights considered
N = 10000 # Number of photons that are traced

# Coronal heights
# h_out = 100 # Maximum height of corona
# N_h = 100 # Number of coronal heights considered
# N = 10000 # Number of photons that are traced

α_vals = [0.1, 0.2, 0.3, 0.4]

# ----------------------------------- #
# Check that all required files exist #
# ----------------------------------- #

begin
    check = true
    for a in a_vals
        for medd in m_edd_full
            if medd == 0.0
                fname = "$file_pref_load-ref-frac-$a-0.0-thin"
            else
                fname = "$file_pref_load-ref-frac-$a-$medd-SS"
            end
            if ! isfile("$load_dir/$fname")
                check = false
                println("File $fname not present in provided directory")
            end 
        end
    end
    if check
        println("All required files present")
    else
        prinln("One or more required files were not detected, ensure presence before proceeding")
    end
end

# ---------------------------- #
# Max R as a function of m_edd
# ---------------------------- #



begin
    p = plot(xlabel="m_edd", ylabel="R", title="Maximum R as a function of m_edd", legend=:bottomright, legendfontsize=7)
    for a_sel in a_vals
        m = KerrMetric(1.0, a_sel)
        heights = create_heights(m, h_out, N_h)
        heights_fine = create_heights(m, h_out, 500)        
        
        R_max = zeros(size(m_edd_full))
        for (n,medd_sel) in enumerate(m_edd_full)
            # Load saved fractions then calculate R
            if medd_sel == 0.0
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-0.0-thin")
            else
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-$medd_sel-SS")
            end
            m_R, m_R_h = calc_max_R(d, heights, heights_fine; ISCO=Gradus.isco(m))
            R_max[n] = m_R                        
        end
        plot!(m_edd_full, R_max, label="a = $a_sel", ls=:dash, markershape=:circle, ms=2.5)    
    end
    yaxis!(ylims=(1,50), :log10, yticks=([1, 2, 5, 10, 20,  50], [1, 2, 5, 10, 20,  50]))
    !isfile(joinpath(save_dir, "R_vs_medd.pdf")) || overwrite ? savefig(joinpath(save_dir, "R_vs_medd.pdf")) : println("Figure already present, not saving")
    display(p)
end

# ------------------------------ #
# Thin disk reflection fractions #
# ------------------------------ #

begin
pl = plot(grid=true, minorgrid=true)
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
!isfile(joinpath(save_dir, "R_ThinDisk.pdf")) || overwrite ? savefig(joinpath(save_dir, "R_ThinDisk.pdf")) : println("Figure already present, not saving")
display(pl)
end

# -------------------------- #
# Source height at maximum R #
# -------------------------- #

begin
    
    Ignore_below_ISCO = false

    max_R = zeros(size(a_vals)) # Stores max R for given spin
    max_R_h = zeros(size(a_vals)) # Stores height at which max R is reached
    ISCO_vals = zeros(size(a_vals))
    EH_vals = zeros(size(a_vals))

    plot()
    for m_edd_sel in m_edd_full[1:1]
        for (n, a) in enumerate(a_vals)    
            
            # Define m to create heights
            m = KerrMetric(1.0, a)
            heights = create_heights(m, h_out, N_h)
            heights_fine = create_heights(m, h_out, 500)

            # Load saved fractions then calculate R
            if m_edd_sel == 0.0
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a-$m_edd_sel-thin")
            else
                d = loaddata(load_dir, "$file_pref_load-ref-frac-$a-$m_edd_sel-SS")
            end
            
            m_R, m_R_h = calc_max_R(d, heights, heights_fine)
            
            max_R[n] = m_R
            max_R_h[n] = m_R_h
        end
        scatter!(a_vals, max_R_h, label="m_edd = $m_edd_sel", shape=:circle, ms=2.5, alpha=1)
    end

    # Plot ISCO and EH
    n = 100
    ISCO_vals = zeros(n)
    EH_vals = zeros(n)
    a_vals_ISCO = LinRange(0, 0.998, n)
    for (n,a) in enumerate(a_vals_ISCO)
        m_temp = KerrMetric(1.0, a)
        ISCO = Gradus.isco(m_temp)
        ISCO_vals[n] = ISCO
        EH_vals[n] = Gradus.inner_radius(m_temp)
    end

    
    # scatter!(a_vals, max_R_h, label="Thin Disk")
    plot!(a_vals_ISCO, ISCO_vals, label="ISCO", ls=:dash, c=:grey)
    plot!(a_vals_ISCO, EH_vals, label="EH", ls=:dash, c=:black)

    xaxis!("Spin")
    yaxis!("Source Height (r_g)", ylims=(1,12), :log10)
    title!("Source height at Maximum R")
    !isfile(joinpath(save_dir, "Rmax_height.pdf")) || overwrite ? savefig(joinpath(save_dir, "Rmax_height.pdf")) : println("Figure already present, not saving")
end


function hmax(a)
    # dauser+2014 fit for max height
    (1.89 * a^2 - 10.86*a + 10.07) * (1 + (9.41e-4 / log10(a)))
end

# Overplot fit from dauser+2014 #
# a_vals_fine = LinRange(0,1, 100)
# hmax_pred = hmax.(a_vals_fine)
# plot!(a_vals_fine, hmax_pred)

# -------------------------- #
# R for different m_edd
# -------------------------- #

begin
    plot()
    for (a, a_sel) in enumerate([0.0, 0.5, 0.998])
        m = KerrMetric(1.0, a_sel)
        for (n, m_sel) in enumerate([0.1, 0.4, 0.7, 1.0])
            heights = create_heights(m, h_out, N_h)
            fname = "$file_pref_load-ref-frac-$a_sel-$m_sel-SS"
            d = loaddata(load_dir, fname)
            plot!(heights, d["above_isco"] ./ d["missed"], label="", ls=linestyles[n], c=colors[a])
            if a == 1
                plot!([], [], ls=linestyles[n], c=:black, label="\$\\dot{m}\$ = $m_sel")
            end
        end
        plot!([], [], ls=:solid, c=colors[a], label="a = $a_sel")
    end
    xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
    yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 25, 50], [1, 2, 5, 10, 25, 50]))
    title!("Reflection Fractions for varying \$\\dot{m}\$")
    !isfile(joinpath(save_dir, "R-vs-h_SS.pdf")) || overwrite ? savefig(joinpath(save_dir, "R-vs-h_SS.pdf")) : println("Figure already present, not saving")
end

# --------------------------- #
# Include Within ISCO photons #
# --------------------------- #

begin
    p = plot()
    plot!([], [], ls=linestyles[1], c=:black, label="With ISCO")
    plot!([], [], ls=linestyles[2], c=:black, label="Without ISCO")
    c_count = 1
    for (a, a_sel) in enumerate([0.0, 0.998])
        m = KerrMetric(1.0, a_sel)
        for (n, m_sel) in enumerate([0.1, 1.0])
            heights = create_heights(m, h_out, N_h)
            fname = "$file_pref_load-ref-frac-$a_sel-$m_sel-SS"
            d = loaddata(load_dir, fname)
            plot!(heights, (d["above_isco"] .+ d["below_isco"]) ./ d["missed"], c=colors[c_count], ls = linestyles[1], label="", markershape=:circle, ms=1)
            plot!(heights, d["above_isco"] ./ d["missed"], label="", c=colors[c_count], ls =linestyles[2], markershape=:circle, ms=1)
            plot!([], [], ls=:solid, c=colors[c_count], label="a = $a_sel, \$\\dot{m}\$ = $m_sel")
            c_count += 1
            if a == 1
                # plot!([], [], ls=linestyles[n], c=:black, label="\$\\dot{m}\$ = $m_sel")
            end
        end
    end
    xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
    yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 25, 50], [1, 2, 5, 10, 25, 50]))
    title!("Reflection Fractions including ISCO photons")
    !isfile(joinpath(save_dir, "R-vs-h-SS-WithISCO.pdf")) || overwrite ? savefig(joinpath(save_dir, "R-vs-h-SS-WithISCO.pdf")) : println("Figure already present, not saving")
    display(p)
end

# -------------------- #
# Funnel Like Geometry #
# -------------------- #

begin
    plot()
    for (a, a_sel) in enumerate([0.0, 0.5, 0.998])
        m = KerrMetric(1.0, a_sel)
        for (n, alpha) in enumerate(α_vals)
            heights = create_heights(m, h_out, N_h)
            fname = "$file_pref-ref-frac-$a_sel-$alpha-funnel"
            d = loaddata(load_dir, fname)
            plot!(heights, d["above_isco"] ./ d["missed"], label="", ls=linestyles[n], c=colors[a])
            if a == 1
                plot!([], [], ls=linestyles[n], c=:black, label="α = $alpha")
            end
        end
        plot!([], [], ls=:solid, c=colors[a], label="a = $a_sel")
    end
xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 25, 50], [1, 2, 5, 10, 25, 50]))
title!("Reflection Fractions for varying α")
savefig("$save_dir/R-vs-h-funnel.pdf")
end


# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #