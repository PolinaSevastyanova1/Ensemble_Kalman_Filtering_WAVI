using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors

ENSEMBLE_NAME = "ensemble_onerun"

@load "$(ENSEMBLE_NAME)/output/eki.jld2" eki param_dict prior

weert, glen, bump, melt, pct = eachrow(eki.u[1].stored_data)

params = [
(weert, "Sliding Law Initialization", "weertman prefactor", "plots/one_hist_weert.png", (0.4,1.7), :grey),
(glen, "Glen's Flow Law Initialization", "glen prefactor", "plots/one_hist_glen.png", (0.4,1.7), :green),
(bump, "Anthropogenic Bump Initialization", "bump prefactor", "plots/one_hist_bump.png", (-100,500), :red),
(melt, "Melt Rate Prefactor Initialization", "melt rate prefactor", "plots/one_hist_melt.png", (0.5,15.5), :blue),
(pct, "Underlying Trend Initialization", "per century trend parameter", "plots/one_hist_pct.png", (-200,400), :cyan)
]

for (param, title, xlabel, filename, xlims, color) in params

	histogram(param, bins=20, xlims=xlims, xlabel=xlabel, ylabel="Frequency", dpi=300, color=color, title=title)
	savefig(filename)
	println("Plot saved in '$filename'")
end

