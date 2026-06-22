using Gradus, Plots
include(joinpath(@__DIR__, "StandardFuncs.jl"))
include(joinpath(@__DIR__, "FunnelGeom.jl"))

const _ROOT = joinpath(@__DIR__, "..")

# -------- #
# Settings #
# -------- #

save_dir = joinpath(_ROOT, "data", "results", "pub", "stash") # Set where to write data to
load_dir = joinpath(_ROOT, "data", "results", "pub", "stash") # Set where to load data from
file_pref = "run_2" # Add a prefix to the filename being saved/loaded
overwrite = false

# Initialize grid

# Spin values for which results are calculated / loaded from
a_vals = [0.0, 0.1, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.900, 0.95, 0.990, 0.998]
# a_vals = [0.0]
# a_vals = [0.0, 0.3, 0.5, 0.7]

# Eddington ratios for which results are calculated / loaded from
m_edd_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 
# m_edd_vals = []

α_vals = [0.1, 0.2, 0.3, 0.4]

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

    if !isfile("$save_dir/$fname") || overwrite
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

        if !isfile("$save_dir/$fname") || overwrite
            d = ShakuraSunyaev(m; eddington_ratio=m_edd) ∘ ThinDisc(0.0, Inf)
            heights, geods = calc_geods(m, d; N,h_out, N_h)
            cf = count_fractions(geods, Gradus.isco(m))
            stashdata(save_dir, fname; cf...)
        else
            println("Data already present, skipping S&S disk")
        end
    end
    
    for alpha in α_vals
        @info "α = $alpha"
        fname = "$file_pref-ref-frac-$a-$alpha-funnel"
        
        if 1 == 1
            d = funnel_disk(alpha, m;) ∘ ThinDisc(0.0, Inf)
            heights, geods = calc_geods(m, d; N, h_out, N_h)
            cf = count_fractions(geods, Gradus.isco(m))
            stashdata(save_dir, fname; cf...)
        else
            println("Data already present, skipping S&S disk")
        end
    end
end