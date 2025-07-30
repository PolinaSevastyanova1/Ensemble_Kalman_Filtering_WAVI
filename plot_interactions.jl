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
members = 10
params = 5
success = false 
coordinates_all = Array{Float64}(undef, params, members)
mean_truth = Float64(mean(y))
std_truth = Float64(std(y))
max_truth = mean_truth + 3*std_truth
min_truth = mean_truth - 3*std_truth

#defining dot colours for too melted, not melted enough, and just right
burgundy_red = RGB(0.5, 0.0, 0.125)
navy_blue = RGB(0.0, 0.0, 0.5)       
bright_green = RGB(0.0, 1.0, 0.0)

#params stored in 6*10 array, of parameters and members, for each iteration
weert, glen, bump, melt, pct = [Array{Float64}(undef, members) for _ in 1:5]

#initialising 3D plot with all members at each iteration, different combinations
plotlyjs()
plot1 = plot(title="3D Parameter Plot 1", seriestype=:scatter, xlabel="melt", ylabel="glen", zlabel="weert", legend=false, size=(400,400))
plot2 = plot(title="3D Parameter Plot 2", seriestype=:scatter, xlabel="melt", ylabel="pct", zlabel="bump", legend=false, size=(400,400))
plot3 = plot(title="3D Parameter Plot 3", seriestype=:scatter, xlabel="glen", ylabel="pct", zlabel="bump", legend=false, size=(400,400))
plot4 = plot(title="3D Parameter Plot 4", seriestype=:scatter, xlabel="weert", ylabel="pct", zlabel="bump", legend=false, size=(400,400))

# iterating over all iterations
for i = 1:iterations
	#extracting param rows and member columns for each iterations
	u_iter = eki.u[i].stored_data
	weert[:], glen[:], bump[:], melt[:],pct[:] = eachrow(u_iter)

	#defining dot colour
	ice_masses = eki.g[i].stored_data
	mean_ice_masses = mean(ice_masses, dims=1)

	colours = [if mean_ice_mass > max_truth 
			navy_blue
		elseif mean_ice_mass < min_truth
			burgundy_red
		else 
			bright_green
		end for mean_ice_mass in mean_ice_masses]

		#plotting in 3D the coordinates of 3 chosen iterations
	plot!(plot1, melt, glen, weert, seriestype=:scatter, xlabel="melt", ylabel="glen", zlabel="weert", title="3D Parameter Plot", marker=:circle,color=colours[:], markersize=2, markerstrokewidth=0)
	plot!(plot1, melt, glen, weert, seriestype=:line, xlabel="melt", ylabel="glen", zlabel="weert", title="3D Parameter Plot", lw=1,color=:grey, alpha=0.8)

	plot!(plot2, melt, pct, bump, seriestype=:scatter, xlabel="melt", ylabel="pct", zlabel="bump", title="3D Parameter Plot", marker=:circle, color=colours[:],markersize=2, markerstrokewidth=0)
	plot!(plot2, melt, pct, bump, seriestype=:line, xlabel="melt", ylabel="pct", zlabel="bump", title="3D Parameter Plot", lw=1,color=:grey,alpha=0.8)

	plot!(plot3, glen, pct, bump, seriestype=:scatter, xlabel="glen", ylabel="pct", zlabel="bump", title="3D Parameter Plot", marker=:circle, color=colours[:],markersize=2, markerstrokewidth=0)
	plot!(plot3, glen, pct, bump, seriestype=:line, xlabel="glen", ylabel="pct", zlabel="bump", title="3D Parameter Plot", lw=1,color=:grey, alpha=0.8)

	plot!(plot4, weert, pct, bump, seriestype=:scatter, xlabel="weert", ylabel="pct", zlabel="bump", title="3D Parameter Plot", marker=:circle, color=colours[:],markersize=2, markerstrokewidth=0)
	plot!(plot4, weert, pct, bump, seriestype=:line, xlabel="weert", ylabel="pct", zlabel="bump", title="3D Parameter Plot", lw=1,color=:grey, alpha=0.8)

end

savefig(plot1, "plots/3DPlot1.html")
savefig(plot2, "plots/3DPlot2.html")
savefig(plot3, "plots/3DPlot3.html")
savefig(plot4, "plots/3DPlot4.html")
println("Plots savedin plots")

