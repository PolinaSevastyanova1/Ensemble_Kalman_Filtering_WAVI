using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors
using Statistics

FOLDER = "ensemble_onerun"
@load "$(FOLDER)/output/truth.jld2" y
obs = y
# add 4 truth values at defined years
scatter_times = [285, 290, 295, 300]

function fit_line(x,y)

	n=length(x)
	x_mean = mean(x)
	y_mean = mean(y)
	slope = sum((x .- x_mean) .* (y .- y_mean)) / sum((x .- x_mean).^2)
	intercept = y_mean - slope * x_mean
	return slope, intercept
end

p = plot(scatter_times, obs, lw=1, color=:red)

scatter!(scatter_times, obs, color=:red, marker=:circle, title="Generated Ice Mass Values", xlabel="Time (years)", ylabel="Ice Mass (10^12 tons)", lw=2, dpi=300, xlims=(250,350), ylims=(10,16))

slope, intercept = fit_line(scatter_times, obs)
best_fit_x = [250,300,350]
best_fit_obs = slope .* best_fit_x .+ intercept
plot!(best_fit_x, best_fit_obs, lw=1, linestyle=:dash, color=:red, legend=:false)

# save plot
savefig(p, "plots/truth_ice_mass_vs_time.png")  # Save the plot to a file
println("Plot saved as 'truth_ice_mass_vs_time.png'")


