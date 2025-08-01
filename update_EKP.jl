include("shared.jl") # packages

function main()

    output_dir = ARGS[1]
    eki_path = joinpath(output_dir, ARGS[2])

    # Parameters
    iteration = parse(Int64, ARGS[3])

    @info "Updating EKP parameters in $(iteration) to $(output_dir)"

    # load current state 
    @load eki_path eki param_dict prior
    N_ensemble = eki.N_ens
    #dim_output = size(eki.observation_series.observations)[1] # size(eki.obs_mean)[1]
    dim_output = 4
    # load data from the ensemble
    G_ens = zeros(dim_output, N_ensemble)
    for member in 1:N_ensemble
        member_path = path_to_ensemble_member(output_dir, iteration, member)
        @load joinpath(member_path, "output.jld2") model_output
        G_ens[:, member] = model_output
    end

    # perform the update    
    EKP.update_ensemble!(eki, G_ens)

    # save the parameter ensemble and EKP
    save_parameter_ensemble(
        get_u_final(eki), # constraints applied when saving
        prior,
        param_dict,
        output_dir,
        "parameters",
        iteration + 1, #save for next iteration
    )

    #save new state
    @save eki_path eki param_dict prior

end

main()
