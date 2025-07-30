using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors
using Statistics
using Distributions
using LaTeXStrings

default(
  fontfamily   = "Computer Modern",  # match LaTeX serif
  titlefont    = font(14,  "Computer Modern", :bold),
  guidefont    = font(14,  "Computer Modern"),
  tickfont     = font(12,  "Computer Modern"),
  legendfont   = font(12,  "Computer Modern"),
  framestyle   = :box,
  grid         = false,
  dpi          = 500,
)


params = [
  (:weert, "Sliding Law Evolution",        "Weertman prefactor",          "hist_hpc_pts/all/comb_weert.png", (0.4,1.7)),
  (:glen,  "Glen's Flow Law Evolution",    "Glen prefactor",              "hist_hpc_pts/all/comb_glen.png",  (0.4,1.7)),
  (:bump,  "Anthropogenic Bump Evolution", "Bump Amplitude (m)",              "hist_hpc_pts/all/comb_bump.png",  (-100,500)),
  (:melt,  "Melt Rate Prefactor Evolution","Melt rate prefactor",         "hist_hpc_pts/all/comb_melt.png",  (0.5,15.5)),
  (:pct,   "Underlying Trend Evolution",   "Anthropogenic trend (m/century)","hist_hpc_pts/all/comb_pct.png",(-300,300))
]

const N_ITER = 20
all_failed  = Dict{Symbol,Array{Float64,2}}()
all_success = Dict{Symbol,Array{Float64,2}}()
all_prior   = Dict{Symbol,Vector{Float64}}()
for (sym,_,_,_,_) in params
  all_failed[sym]  = Array{Float64}(undef, N_ITER, 0)
  all_success[sym] = Array{Float64}(undef, N_ITER, 0)
  all_prior[sym]   = Float64[]
end


FAIL_COLOR    = "#4A90E2"  # light blue
SUCCESS_COLOR = "#800020"  # burgundy red
const COL0  = "#4A90E2"
@load "ensemble/output/eki.jld2" eki param_dict prior
@load "ensemble/output/truth.jld2" y
global first_failed    = true
global first_success   = true
global first_reference = true
global obs = eki.observation_series.observations[1].samples
μ_truth  = obs[1]
mean_truth = 12.52
tol           = 1.0

# initialising arrays
iterations = length(eki.u)-1
members = 20
error_values = Array{Float64}(undef, iterations, members)
mean_truth = Float64(mean(y))
std_truth = Float64(std(y))
# iterating over all iterations
for i = 1:iterations
	g_matrix = eki.g[i].stored_data
	#mean_array = mean(g_matrix, dims=1)
	#error_values[i,:] = mean_array
	error_values[i, :] = g_matrix[end, :] 
end
g_final = eki.g[end].stored_data
mask_ok = [ all(abs.(g_final[:, m] .- μ_truth) .<= tol) for m in 1:members ]
ok_ix   = findall(mask_ok)
fail_ix = findall(!,    mask_ok)
plot(xlabel="Iteration", ylabel="Ice Mass " * L"( \times 10^{12}\,\mathrm{t})", dpi=500, legend=:topright)

for member in 1:members
	lbl = first_failed ? "Member (failed)" : false
	plot!(1:iterations, error_values[:,member], label=lbl, lw=1, color=FAIL_COLOR, alpha = 0.3)
	global first_failed = false
end


for member in ok_ix
    lbl = first_success ? "Member (successful)" : false
    plot!(
        1:iterations,
        error_values[:, member],
        lw    = 1,
        color = SUCCESS_COLOR,
        alpha = 0.8,
        label = lbl
    )
    global first_success = false
end

mean_truth = 12.52
# 3) reference lines, but only label once
if first_reference
        hline!([mean_truth],           linestyle=:dash, color=:black, label = "Mean ice volume")
        hline!([mean_truth + tol],     linestyle=:dash, color=:grey,  label="Tolerance ranges")
        hline!([mean_truth - tol],     linestyle=:dash, color=:grey,  label=false)
        global first_reference = false
else
        hline!([mean_truth],           linestyle=:dash, color=:black, label=false)
        hline!([mean_truth + tol],     linestyle=:dash, color=:grey,  label=false)
        hline!([mean_truth - tol],     linestyle=:dash, color=:grey,  label=false)
end

savefig("plots/ice_volume_mean_iterations.png")
println("Plot saved in plots/ice_volume_mean_iterations.png")


