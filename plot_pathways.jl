# using JLD2
# using EnsembleKalmanProcesses
# using NetCDF
# using Plots
# using Printf
# using Colors


# FOLDER = "ensemble"
# members = 20
# path = "$(FOLDER)/output/eki.jld2"
# @load path eki param_dict prior
# obs = eki.observation_series.observations[1].samples
# @load "$(FOLDER)/output/truth.jld2" y
# ITERATE_TO = length(eki.u)-1

# # initilising plot
# p = plot(title="Ice Mass vs Time", xlabel="Time (years)", ylabel="Ice Mass (10^12 tons)", lw=2, dpi=300)

# # light blue to dark blue gradient
# colors = [RGB(153/255, 153/255, 255/255), RGB(102/255, 102/255, 255/255), RGB(0, 0, 204/255), RGB(0, 0, 153/255), RGB(0, 0, 65/255)]

# #iterating over iterations
# for j in 0:ITERATE_TO
#     iteration_path = @sprintf("%s/output/iteration_%03d", FOLDER, j)
#     color = colors[1]
#     # color = colors[trunc(Int, j/2)+1]
    
#     if j==ITERATE_TO
#         color = :red 
#     end  # last iteration is red

#     # iterating over members
#     for i in 1:members
#         member_path = @sprintf("%s/member_%03d/outfile.nc", iteration_path, i)
#         nc_data = NetCDF.open(member_path)
#         h_data = nc_data["h"]
#         time = nc_data["TIME"][:]

#         # ice mass for memeber
#         ice_mass = zeros(Float64, size(h_data, 3))
#         for t = 1:size(h_data, 3)
#             ice_mass[t] = sum(h_data[:, :, t]) * 2.0e3 * 2.0e3 / 1e12
#         end

#         # add member ice mass to plot
#         plot!(p, time, ice_mass, label="", color=color)  # Same color for all members
#     end
# end

# # add 4 truth values at defined years
# scatter_times = [285, 290, 295, 300]
# scatter!(p, scatter_times, obs, color=:red, marker=:circle, label="Truth")


# # save plot
# savefig(p, "plots/ice_mass_vs_time.png")  # Save the plot to a file
# println("Plot saved as 'ice_mass_vs_time.png'")

using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Statistics
using LaTeXStrings
using Printf
using Measures


# ──────────────────────────────────────────────────────────────────────────────
# 1) GLOBAL STYLING
# ──────────────────────────────────────────────────────────────────────────────
default(
  fontfamily   = "Computer Modern",
  titlefont    = font(14, "Computer Modern", :bold),
  guidefont    = font(14, "Computer Modern"),
  tickfont     = font(12, "Computer Modern"),
  legendfont   = font(12, "Computer Modern"),
  framestyle   = :box,
  grid         = false,
  dpi          = 500,
)

# ──────────────────────────────────────────────────────────────────────────────
# 2) DATA FOLDERS & CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────
folders    = ["EKFpt13"]
NMEM  = 20
OUTPN = "plots/pathways.png"
FAIL_COLOR    = "#4A90E2"  # light blue
SUCCESS_COLOR = "#800020"  # burgundy red
const COL0  = "#4A90E2"

# ──────────────────────────────────────────────────────────────────────────────
# 3) INITIALIZE PLOT
# ──────────────────────────────────────────────────────────────────────────────
p = plot(
  title  = "",
  xlabel = "Time (years)",
  ylabel = "Ice Mass " * L" (\times 10^{12}\,\mathrm{t})",
  legend = :topleft,
    legendfontsize  = 9,
  legendframe     = :box,
  legendbg        = :white,
  legendborder    = :black,
  dpi    = 500,
  margin = 3mm
)

# ──────────────────────────────────────────────────────────────────────────────
# 4) ACCUMULATE FINAL-ITERATION PATHWAYS
# ──────────────────────────────────────────────────────────────────────────────
t0_list = Vector{Vector{Float64}}()
im0_list = Vector{Vector{Float64}}()
tseries_list = Vector{Vector{Float64}}()
imass_list   = Vector{Vector{Float64}}()
success_list = Bool[]
wide_success_list = Bool[]

for folder in folders
    @load "ensemble/output/eki.jld2" eki
    @load "ensemble/output/truth.jld2" y

    global obs = eki.observation_series.observations[1].samples
    # tolerance criteria
    # μ_truth = mean(y)
    # tol     = 0.5 * std(y)

    μ_truth  = obs[1]
    tol     = 1.0
    wide_tol = 2.0
    obs_times = [285, 290, 295, 300]

    # final iteration index
    #JFINAL = 9
    JFINAL = folder in ("EKFpt16","EKFpt17") ? 7 : 19
    for m in 1:NMEM
        path0 = @sprintf("ensemble/output/iteration_%03d", 0)
        file0 = @sprintf("%s/member_%03d/outfile.nc", path0, m)
        if isfile(file0)
            NetCDF.open(file0) do ds
                tseries = ds["TIME"][:]
                hdat     = ds["h"]
                imass0   = [sum(hdat[:,:,t]) * 2e3 * 2e3 / 1e12 for t in 1:size(hdat,3)]
                push!(t0_list, tseries)
                push!(im0_list, imass0)
            end
        else
            @warn "Missing file: $file0"
        end

        # final iteration
        iter_path = @sprintf("ensemble/output/iteration_%03d", JFINAL)
        ncfile    = @sprintf("%s/member_%03d/outfile.nc", iter_path, m)

        if !isfile(ncfile)
            @warn "Missing file: $ncfile"
            continue
        end

        NetCDF.open(ncfile) do ds
            hdat    = ds["h"]
            tseries = ds["TIME"][:]
            imass   = [ sum(hdat[:,:,t]) * 2e3 * 2e3 / 1e12
                        for t in 1:size(hdat,3) ]

            push!(tseries_list, tseries)
            push!(imass_list,   imass)
            #push!(success_list, abs(imass[end] - μ_truth) <= tol)
            #push!(success_list, abs(imass[285, 290, 295, 300] - μ_truth) <= tol)
            indices = [findfirst(==(t), tseries) for t in obs_times]
            model_vals = imass[indices]
            is_success = all(abs.(model_vals .- μ_truth) .<= tol)
            is_wide_success = all(abs.(model_vals .- μ_truth) .<= wide_tol)
            push!(success_list, is_success)
            push!(wide_success_list, is_wide_success)
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# 5) PLOT FAILED (light blue) THEN SUCCESSFUL (burgundy red)
# ──────────────────────────────────────────────────────────────────────────────
for i in eachindex(im0_list)
    plot!(p, t0_list[i], im0_list[i];
          color = COL0, alpha = 0.6, lw = 1.5, label = false)
end

for i in eachindex(imass_list)
    if !success_list[i]
        plot!(p, tseries_list[i], imass_list[i];
              color = :indianred, alpha = 0.5, lw = 1.5, label = false)
    end
end

for i in eachindex(imass_list)
    if success_list[i]
        plot!(p, tseries_list[i], imass_list[i];
              color = SUCCESS_COLOR, alpha = 0.9, lw = 1.5, label =false)
    end
end

println("Found $(length(success_list)) members, of which $(count(success_list)) were successful.")
println("Found $(length(wide_success_list)) members, of which $(count(wide_success_list)) were successful within 2*err.")
# ──────────────────────────────────────────────────────────────────────────────
# add 4 truth values at defined years
scatter_times = [285, 290, 295, 300]
scatter!(
  p,
  scatter_times,
  obs;
  yerror            = fill(1.0, length(obs)),     # ±1 error bars
  color             = "#800020",                  # same burgundy as your successes
  alpha             = 1.0,
  marker            = :circle,
  ms                = 4,                          # slightly larger markers
  markerstrokecolor = :black,                     # crisp black outline
  markerstrokewidth = 1.3,
  lw                = 1,                          # error‐bar line‐width
  label             = "Truth"                    # LaTeXString to match fonts
)

# dummy series just to drive the legend
plot!(
  p,
  [NaN], [NaN];
  color = COL0,
  lw    = 2,
  label = "Prior"
)
plot!(
  p,
  [NaN], [NaN];
  color = :indianred,
  lw    = 2,
  label = "Posterior (failed)"
)
plot!(
  p,
  [NaN], [NaN];
  color = SUCCESS_COLOR,
  lw    = 2,
  label = "Posterior (successful)"
)

plot!(
  p;
  xticks = (0:50:300, string.(1700:50:2000)),
  xlim   = (0, 300),
  ylim = (0,35)
)


# ──────────────────────────────────────────────────────────────────────────────
# 6) SAVE
# ──────────────────────────────────────────────────────────────────────────────
savefig(p, OUTPN)
println("Saved combined pathways plot to '$(OUTPN)'")
