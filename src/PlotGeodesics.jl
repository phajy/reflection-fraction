using Gradus
using Plots
include("StandardFuncs.jl")
include("CustomShakuraSunyaev.jl")
using RecipesBase

# --------------------------------------------------------------------------------
# !!!!!!
# To plot the traced rays, the "save" option in tracegeodesics needs to be 
# set to true in StandardFuncs.jl
# The solver also needs to be commented out so that it uses the default one.
# (this is a bit inconvenient but I don't want to take the time to code this up correctly)
# !!!!!!
# --------------------------------------------------------------------------------

# Calculate Geodesics
begin
a = 0.0
h_out = 100
N_h = 2
thrshld = 0.1
m_edd = 0.3
N=1000

m = KerrMetric(1.0, a)

d = ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld)
heights, geods_custom = calc_geods(m, d; N,h_out, N_h)

d = ThinDisc(0.0, 500.0)
heights, geods_thin = calc_geods(m, d; N,h_out, N_h)

print(":)") # To ensure that it doesn't print out all the returned geodesics :)
end

# ---------------------------------------------------#
# Plot all traced traced paths that reach infinity   #
# ---------------------------------------------------#

height_idx = 2

# Thick Disk
sol = geods_custom[height_idx].gps
itr = sol
plot(legend=:topleft, ylabel="Height (rg)", xlabel="- radius (rg)")
coord_max = nothing
for s in itr.u
    if Gradus.get_status_code(s.prob.p) == StatusCodes.NoStatus || Gradus.get_status_code(s.prob.p) == StatusCodes.OutOfDomain
        temp = Gradus._extract_path(s, 1000; t_span=10000.0)
        coord = (temp[1], temp[3])
        if isnothing(coord_max)
             plot!(coord[1], coord[2], ls=:solid, c=:red, alpha=0.5, label="")
        end
        coord_max = coord
        plot!(coord[1], coord[2], ls=:dot, label="", c=:red)
    end
end
# coord_max
# plot!(coord_max[1], coord_max[2], c=:blue2, ls=:solid)

sol = geods_thin[height_idx].gps
itr = sol
coord_max = nothing
for s in itr.u
    if Gradus.get_status_code(s.prob.p) == StatusCodes.NoStatus || Gradus.get_status_code(s.prob.p) == StatusCodes.OutOfDomain
        temp = Gradus._extract_path(s, 1000; t_span=10000.0)
        coord = (temp[1], temp[3])
        plot!(coord[1], coord[2], ls=:dot, label="", c=:blue)
        if isnothing(coord_max)
             plot!(coord[1], coord[2], ls=:solid, c=:blue, alpha=1,  label="")
        end
        coord_max = coord
    end
end
r = LinRange(0.1, 600, 1000)
cs_ss = Gradus.cross_section.(ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld), r)
plot!(-r, cs_ss, c=:black,  label="Disk Cross Section")

xaxis!(xlims=(-550, 0))
yaxis!(ylims=(-10, 200))

# Zoom in
# xaxis!(xlims=(-20, 0))
# yaxis!(ylims=(-10, 100))

title!("Traced photons that reach infinity")
# savefig("data/results/custom_shakura_sunyaev/Debugging/ThickDiskPhotonTrajINf_Thr$thrshld.pdf")


# ---------------------------------------------------#
# Plot paths that intersect the disk
# ---------------------------------------------------#

height_idx = 2

# Thick Disk
sol = geods_custom[height_idx].gps
itr = sol
plot(legend=:topleft, ylabel="Height (rg)", xlabel="- radius (rg)")
coord_max = nothing
for s in itr.u
    if Gradus.get_status_code(s.prob.p) == StatusCodes.IntersectedWithGeometry
        temp = Gradus._extract_path(s, 1000; t_span=10000.0)
        coord = (temp[1], temp[3])
        if isnothing(coord_max)
             plot!(coord[1], coord[2], ls=:solid, c=:blue2, alpha=0.5, label="")
        end
        coord_max = coord
        plot!(coord[1], coord[2], ls=:dot, label="")
    end
end
# coord_max
plot!(coord_max[1], coord_max[2], c=:blue2, ls=:solid, label="")

sol = geods_thin[height_idx].gps
itr = sol
coord_max = nothing
for s in itr.u
    if Gradus.get_status_code(s.prob.p) == StatusCodes.IntersectedWithGeometry
        temp = Gradus._extract_path(s, 1000; t_span=10000.0)
        coord = (temp[1], temp[3])
        plot!(coord[1], coord[2], ls=:dot, label="", alpha=.5)
        if isnothing(coord_max)
             plot!(coord[1], coord[2], ls=:solid, c=:crimson, alpha=.5, label="")
        end
        coord_max = coord
    end
end
# coord_max
plot!(coord_max[1], coord_max[2], c=:crimson, ls=:solid, label="")

# hline!([0.0], color=:black)
# plot_horizon!(m)
xaxis!(xlims=(-550, 0))
yaxis!(ylims=(-10, 130))

# Zoom In
# xaxis!(xlims=(-20, 0))
# yaxis!(ylims=(-10, 100))

r = LinRange(0.1, 600, 1000)
cs_ss = Gradus.cross_section.(ShakuraSunyaev_custom(m; eddington_ratio=m_edd, threshold=thrshld), r)
plot!(-r, cs_ss, c=:black,  label="Disk Cross Section")

hline!([0], c=:black, ls=:dash, label="Thin Disk")

# To add labels for plot of both disc and inf photons
# plot!([], [], c=:red, ls=:dot, label="Infinity")
# plot!([], [], c=:blue, ls=:dot, label="Disc")

title!("Traced photons that hit the disc")
savefig("data/results/custom_shakura_sunyaev/Debugging/PhotonTrajectory_Disc_Thr$thrshld.pdf")
