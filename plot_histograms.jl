using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors

ITER_NO = 30 #make this divisible by 5 please
ENSEMBLE_NAME = "ensemble_4"

@load "$(ENSEMBLE_NAME)/output/eki.jld2" eki param_dict prior

function get_iter_params(EKI, iter)
	return EKI.u[iter].stored_data
end

#chunks of 5 iterations, so 50 members in each histogram
chunks_no = Integer(ITER_NO/5)-1
for ch = 0:chunks_no
	i_start = ch*5 +1
	i_end = ((ch+1)*5)
	# initialize empty arrays for storing aggregated
	weertman_params_chunk, glen_params_chunk, bump_params_chunk, melt_params_chunk, pct_params_chunk = [[],[],[],[],[]]

	# iterating over 5 iterations and collecting all parameter values from the chunk
	for i = i_start:i_end
		iter_params = get_iter_params(eki,i)
		weertman_params, glen_params, bump_params,
melt_params, pct_params = eachrow(iter_params)
		append!(weertman_params_chunk, weertman_params)
		append!(glen_params_chunk, glen_params)
		append!(bump_params_chunk, bump_params)
		append!(melt_params_chunk, melt_params)
		append!(pct_params_chunk, pct_params)
	end

	histogram(weertman_params_chunk, bins=25, xlims=(0.4,1.7), color=:grey, title="Sliding Law Parameters, Group $(ch)", xlabel="Prefactor", ylabel="Frequency")
	savefig("histograms/weertman_c_ch$(ch)_hist.png")
	println("Histogram saved in histograms/weertman_c_ch$(ch)_hist.png")

	histogram(glen_params_chunk, bins=25, xlims=(0.4,1.7), color=:green, title="Viscosity Parameters, Group $(ch)", xlabel="Prefactor", ylabel="Frequency")
	savefig("histograms/glen_ch$(ch)_hist.png")
	println("Histogram saved in histograms/glen_ch$(ch)_hist.png")

	histogram(bump_params_chunk, bins=25, xlims=(-100,500), color=:red, title="Anthropogenic Bump Parameters, Group $(ch)", xlabel="Bump Value", ylabel="Frequency")
	savefig("histograms/bump_ch$(ch)_hist.png")
	println("Histogram saved in histograms/bump_ch$(ch)_hist.png")

	histogram(melt_params_chunk, bins=25, xlims=(0.5,15.5), color=:blue, title="Melt Prefactor Parameters, Group $(ch)", xlabel="Prefactor", ylabel="Frequency")
	savefig("histograms/melt_ch$(ch)_hist.png")
	println("Histogram saved in histograms/melt_ch$(ch)_hist.png")

	histogram(pct_params_chunk, bins=25, xlims=(-200,400), title="Per Century Trend Parameters, Group $(ch)", xlabel="Prefactor", ylabel="Frequency")
	savefig("histograms/pct_ch$(ch)_hist.png")
	println("Histogram saved in histograms/pct_ch$(ch)_hist.png")

end


