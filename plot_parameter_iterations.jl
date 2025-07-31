using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors
using Statistics

@load "ensemble/output/eki.jld2" eki param_dict prior

# initialising arrays
iterations = length(eki.u)-1
members = 20
success_mems = []
fail_mems = []
error_values = Array{Float64}(undef, iterations, members)

#params stored in 6*10 array, of parameters and members, for each iteration
weert, glen, bump, melt, pct = [Array{Float64}(undef, iterations, members) for _ in 1:5]

# iterating over all iterations
for i = 1:iterations
	u_iter = eki.u[i].stored_data
	weert[i, :], glen[i,:], bump[i,:], melt[i,:],pct[i,:] = eachrow(u_iter)
end

#param, plot title, ylabel, figure title, 
params = [
(weert, "Sliding Law Evolution", "weertman prefactor", "plots/weert_param_iters.png", (0.4,1.7)),
(glen, "Glen's Flow Law Evolution", "glen prefactor", "plots/glen_param_iters.png", (0.4,1.7)),
(bump, "Anthropogenic Bump Evolution", "bump prefactor", "plots/bump_param_iters.png", (-100,500)),
(melt, "Melt Rate Prefactor Evolution", "melt rate prefactor", "plots/melt_param_iters.png", (0.5,15.5)),
(pct, "Underlying Trend Evolution", "per century trend parameter", "plots/pct_param_iters.png", (-200,400))
]

#and highlighting members 3,5,8,9 which led to the correct values. 
#6 and 10 failed, so highlighting those as well, in dotted lines

colours = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink,:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
for (param, title, ylabel, filename, ylims) in params

	plot(title=title, xlabel="Iteration no.", ylabel=ylabel, dpi=300, legend=:outerright, ylims=ylims)

	for member in 1:members	
		if member in success_mems
			plot!(1:iterations, param[:,member], label="Member $member", lw=2, color=colours[member])
		elseif member in fail_mems
			plot!(1:iterations, param[:,member], label="Member $member", lw=2, linestyle=:dot, color=colours[member])
		else
			plot!(1:iterations, param[:,member], label="Member $member", lw=1, color=colours[member])
		end
	end

	savefig(filename)
	println("Plot saved in '$filename'")
end



