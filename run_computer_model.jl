include("exp_driver.jl")

using TOML, JLD2

# copied from EKP to avoid additional precompile time 
######
function path_to_ensemble_member(load_path, iteration, member; pad_zeros = 3)

    # Get the directory of the iteration
    load_dir = joinpath(load_path, join(["iteration", lpad(iteration, pad_zeros, "0")], "_"))
    subdir_name = join(["member", lpad(member, pad_zeros, "0")], "_")

    return joinpath(load_dir, subdir_name)

end
get_parameter_values(param_dict::Dict, names) = Dict(n => param_dict[n]["value"] for n in names)
######

function main()
    # Paths
    output_dir = ARGS[1]
    data_path = joinpath(output_dir, ARGS[2])
    iteration = parse(Int64, ARGS[3])
    member = parse(Int64, ARGS[4])

    # get parameters
    member_path = path_to_ensemble_member(output_dir, iteration, member)
    param_dict = TOML.parsefile(joinpath(member_path, "parameters.toml"))
    names = ["weertman_c_prefactor", "glen_a_ref_prefactor", "bump_amplitude", "melt_rate_prefactor", "per_century_trend"]
    params = get_parameter_values(param_dict, names)

    # get rng
    #@load data_path rng_model

    # evaluate map with noise to create data
    #model_output = parameter_to_data_map(params, rng = rng_model)
    model_output = parameter_to_data_map(params,member_path)

    output_path = joinpath(member_path, "output.jld2")
    @info "Outputting to $output_path"
    @save output_path model_output

    #save rng
    #@save data_path rng_model

end

main()
