using Gradus, Plots
include("StandardFuncs.jl")

# -------- #
# Settings #
# -------- #

save_dir = "$(pwd())/data/results/pub/stash" # Set where to write data to
load_dir = "$(pwd())/data/results/pub/stash" # Set where to load data from
file_pref = "test_run" # Add a prefix to the filename being saved/loaded
overwrite = false

# Initialize grid

# Spin values for which results are calculated / loaded from
# a_vals = [0.0, 0.3, 0.5, 0.7, 0.900, 0.990, 0.998]
a_vals = [0.0, 0.3]

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
    fname = "$file_pref-frac-$a-0.0-thin"

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