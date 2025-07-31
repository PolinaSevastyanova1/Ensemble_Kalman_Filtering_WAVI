using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Printf
using Colors


FOLDER = "ensemble"

path = "$(FOLDER)/output/eki.jld2"
@load path eki param_dict prior
obs = eki.observation_series.observations[1].samples
@load "$(FOLDER)/output/truth.jld2" y

ITERATE_TO = 9
# initilising plot
p = plot(title="Ice Mass vs Time", xlabel="Time (years)", ylabel="Ice Mass (10^12 tons)", lw=2, dpi=300)
members = 20
# light blue to dark blue gradient
colors = [RGB(153/255, 153/255, 255/255), RGB(102/255, 102/255, 255/255), RGB(0, 0, 204/255), RGB(0, 0, 153/255), RGB(0, 0, 65/255)]

#iterating over iterations
for j in 0:ITERATE_TO
    iteration_path = @sprintf("%s/output/iteration_%03d", FOLDER, j)
    color = colors[1]
    # color = colors[trunc(Int, j/2)+1]
    
    if j==ITERATE_TO
        color = :red 
    end  # last iteration is red

    # iterating over members
    for i in 1:members
        member_path = @sprintf("%s/member_%03d/outfile.nc", iteration_path, i)
        nc_data = NetCDF.open(member_path)
        h_data = nc_data["h"]
        time = nc_data["TIME"][:]

        # ice mass for memeber
        ice_mass = zeros(Float64, size(h_data, 3))
        for t = 1:size(h_data, 3)
            ice_mass[t] = sum(h_data[:, :, t]) * 2.0e3 * 2.0e3 / 1e12
        end

        # add member ice mass to plot
        plot!(p, time, ice_mass, label="", color=color)  # Same color for all members
    end
end

# add 4 truth values at defined years
scatter_times = [285, 290, 295, 300]
scatter!(p, scatter_times, obs, color=:red, marker=:circle, label="Truth")


# save plot
savefig(p, "plots/ice_mass_vs_time.png")  # Save the plot to a file
println("Plot saved as 'ice_mass_vs_time.png'")
