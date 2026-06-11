using GLMakie
using Makie
using Overlay
using Images
include(joinpath(@__DIR__, "src", "StandardFuncs.jl"))

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

# Eddington ratios for which results are calculated / loaded from
m_edd_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 

# Coronal heights
h_out = 100 # Maximum height of corona
N_h = 40 # Number of coronal heights considered
N = 10000 # Number of photons that are traced

file_pref_load = "run_2"

# load Dauser's figure
img = load("dauser_2014.png")
# calibrate the image (either by hand or load existing calibration)
# either
# cal = calibrate_image(img; title = "Click two x, then two y")
# save_calibration("dauser_2014.calibration", cal)
# or
cal = load_calibration("dauser_2014.calibration")

# plot Dauser's figure
fig = Makie.Figure()
ax = Makie.Axis(fig[1, 1])
plot_calibrated_image!(ax, img, cal)

# Thin disk reflection fractions (same data coordinates as the calibrated image)
colors = Makie.wong_colors()
for (n, a) in enumerate(a_vals)
    n in (2, 4) && continue  # skip a = 0.3 and a = 0.7
    m = KerrMetric(1.0, a)
    # Horizon-based heights (same as src/ref-frac.jl); use if stash data used inner_radius, not ISCO
    heights = collect(logrange(Gradus.inner_radius(m) + 1e-2, h_out, N_h))
    fname = "$file_pref_load-ref-frac-$a-0.0-thin"
    d_thin = loaddata(load_dir, fname)
    y = d_thin["above_isco"] ./ d_thin["missed"]
    Makie.lines!(
        ax, heights, y;
        label = "a = $a",
        color = colors[mod1(n, length(colors))],
        linewidth = 2,
    )
end
ax.xlabel = "Source Height (r_g)"
ax.ylabel = "R"
ax.title = "Reflection Fractions"
ax.xticks = ([2, 10, 30, 100], string.([2, 10, 30, 100]))
ax.yticks = ([1, 2, 5, 10, 20], string.([1, 2, 5, 10, 20]))
ax.xgridvisible = true
ax.ygridvisible = true
ax.xminorgridvisible = true
ax.yminorgridvisible = true

Makie.axislegend(ax; position = :rt)
Makie.display(fig)
