using Gradus, Plots
include("StandardFuncs.jl")

# -------- #
# Settings #
# -------- #

save_dir = "$(pwd())/data/results/pub/stash" # Set where to write data to
load_dir = "$(pwd())/data/results/pub/stash" # Set where to load data from
file_pref = "run_2" # Add a prefix to the filename being saved/loaded
overwrite = false

# Initialize grid

# Spin values for which results are calculated / loaded from
a_vals = [0.0, 0.3, 0.5, 0.7, 0.900, 0.990, 0.998]
# a_vals = [0.0, 0.3, 0.5, 0.7]

# Eddington ratios for which results are calculated / loaded from
m_edd_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 
# m_edd_vals = [0.1]

# Coronal heights
h_out = 100 # Maximum height of corona
N_h = 40 # Number of coronal heights considered
N = 10000 # Number of photons that are traced

# ---------------------- #
# Trace and save results #
# ---------------------- #

for a in a_vals
    @info "a = $a"
    m = KerrMetric(1.0, a)
    fname = "$file_pref-ref-frac-$a-0.0-thin"

    if !isfile("$save_dir/$file_pref-frac-$a-0.0-thin") || overwrite
            d_thin = ThinDisc(0.0, Inf)
            heights_thin, geods_thin = calc_geods(m, d_thin; N,h_out, N_h)
            cf = count_fractions(geods_thin, Gradus.isco(m))
            stashdata(save_dir, fname; cf...)
        else
            println("Data already present, skipping thin disk")
    end

    for m_edd in m_edd_vals
        
        @info "m_edd = $m_edd"
        fname = "$file_pref-ref-frac-$a-$m_edd-SS"

        if !isfile("$save_dir/$file_pref-ref-frac-$a-$m_edd-SS") || overwrite
            d = ShakuraSunyaev(m; eddington_ratio=m_edd) ∘ ThinDisc(0.0, Inf)
            heights, geods = calc_geods(m, d; N,h_out, N_h)
            cf = count_fractions(geods, Gradus.isco(m))
            stashdata(save_dir, fname; cf...)
        else
            println("Data already present, skipping S&S disk")
        end
    end
end


# -------------------- #
# Plot from saved data #
# -------------------- #
# file_pref_load = "test_run_loc"

# # Thin disk reflection fractions
# begin
# pl = plot(grid=true, minorgrid=true)
# # colors = [:red, :blue, :green, :cyan, :purple, :orange]
# for (n, a) in enumerate(a_vals)
    
#     m = KerrMetric(1.0, a)
#     heights = create_heights(m, h_out, N_h)

#     fname = "$file_pref_load-ref-frac-$a-0.0-thin"
#     d_thin = loaddata(load_dir, fname)

#     plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="a = $a", ls=:solid) 
# end
# xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
# yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20], [1, 2, 5, 10, 20]))
# title!("Reflection Fractions")
# display(pl)
# end


# begin
# pl = plot(grid=true, minorgrid=true)
# # colors = [:red, :blue, :green, :cyan, :orange, :purple, "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB"]

# colors = [:blue, :black, :maroon, :green, :pink, :purple, :teal, :orange, :red, :silver]
# # Set which accretion rate is being plotted
# m_edd_sel = 1.0

# plot!([], [], label="Thin Disk", ls=:dash, c=:black)
# plot!([], [], label="M_edd = $m_edd_sel", ls=:solid, c=:black)

# for (n, a) in enumerate(a_vals)
    
#     m = KerrMetric(1.0, a)
#     heights = create_heights(m, h_out, N_h)

#     d_custom = loaddata(load_dir, "$file_pref_load-ref-frac-$a-$m_edd_sel-SS")

#     plot!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="a = $a", ls=:solid, c=colors[n]) 

#     d_thin = loaddata(load_dir, "$file_pref_load-ref-frac-$a-0.0-thin")

#     plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="", ls=:dash, c=colors[n]) 
#     # break
# end

# xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
# yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20, 25], [1, 2, 5, 10, 20, 25]))
# title!("Reflection Fractions, \$\\dot{m} = $m_edd_sel\$")
# # savefig("data/results/custom_shakura_sunyaev/Preliminary/R_comparison_medd$m_edd_sel.pdf")
# display(pl)
# end



# begin
# pl = plot(grid=true, minorgrid=true)
# # colors = [:red, :blue, :green, :cyan, :orange, :purple, "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB", "#67BEDB"]

# colors = [:blue, :black, :maroon, :green, :pink, :purple, :teal, :orange, :red, :silver]
# # Set which accretion rate is being plotted
# a_sel = 0.99

# # plot!([], [], label="Thin Disk", ls=:dash, c=:black)
# # plot!([], [], label="M_edd = $m_edd_sel", ls=:solid, c=:black)


# for (n, m_edd) in enumerate(m_edd_vals)
    
#     m = KerrMetric(1.0, a_sel)
#     heights = create_heights(m, h_out, N_h)

#     d_custom = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-$m_edd-SS")

#     plot!(heights, d_custom["above_isco"] ./ d_custom["missed"], label="m_edd = $m_edd", ls=:solid, c=colors[n])     

#     if n == 1
#         d_thin = loaddata(load_dir, "$file_pref_load-ref-frac-$a_sel-0.0-thin")
#         plot!(heights, d_thin["above_isco"] ./ d_thin["missed"], label="m_edd = 0.0", ls=:solid, c=colors[end]) 
#     end
#     # break
# end

# xaxis!("Source Height (r_g)", :log10, xticks=([2, 10, 30, 100], [2, 10, 30, 100]))
# yaxis!("R", :log10, minorgrid=true, yticks=([1, 2, 5, 10, 20, 25], [1, 2, 5, 10, 20, 25]))
# title!("Reflection Fractions, \$ a = $a_sel\$")
# # savefig("data/results/custom_shakura_sunyaev/Preliminary/R_comparison_medd$m_edd_sel.pdf")
# display(pl)
# end