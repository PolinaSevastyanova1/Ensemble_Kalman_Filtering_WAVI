using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors
using Statistics

@load "ensemble/output/eki.jld2" eki param_dict prior
@load "ensemble/output/truth.jld2" y

# initialising arrays
iterations = length(eki.u)-1
members = 20
error_values = Array{Float64}(undef, iterations, members)
mean_truth = Float64(mean(y))
std_truth = Float64(std(y))

# iterating over all iterations
for i = 1:iterations
	g_matrix = eki.g[i].stored_data
	mean_array = mean(g_matrix, dims=1)
	error_values[i,:] = mean_array
end

plot(title="Ensemble member ice volume mean with iteration", xlabel="Iteration no.", ylabel="Ice Volume", dpi=300, legend=:outerright)
colours = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink, :red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]

for member in 1:members
	plot!(1:iterations, error_values[:,member], label="Member $member", lw=1, color=colours[member])
end

hline!([mean_truth], linestyle=:dash, color=:black, label="true ice volume")
hline!([mean_truth] + 0.5*[std_truth], linestyle=:dash, color=:grey, label="true ice volume range")
hline!([mean_truth] - 0.5*[std_truth], linestyle=:dash, color=:grey)

savefig("plots/ice_volume_mean_iterations.png")
println("Plot saved in plots/ice_volume_mean_iterations.png")



