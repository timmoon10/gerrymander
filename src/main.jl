include(joinpath(@__DIR__, "Gerrymander.jl"))

function main(args::Vector{String})

    # Get states from command-line arguments
    state_ids = Set{UInt}()
    state_names = Gerrymander.DataFiles.state_names()
    if isempty(args)
        push!(state_ids, 6)  # Default is California
    elseif length(args) == 1 && args[1] == "all"
        union!(state_ids, keys(state_names))
    else
        for arg in args
            push!(state_ids, parse(UInt, arg))
        end
    end
    state_ids = collect(state_ids)
    sort!(state_ids)

    # Print states
    message = "States: "
    for (i, state_id) in enumerate(state_ids)
        if i > 1
            message *= ", "
        end
        message *= state_names[state_id]
    end
    println(message)

    # Partitioner
    num_partitions::UInt = 6
    personal_stdev::Float64 = 100.0
    max_distance::Float64 = 400.0
    partitioner = Gerrymander.SimulatedAnnealing.Partitioner(
        num_partitions,
        state_ids,
        personal_stdev,
        max_distance,
    )

    # Plotter
    plotter = Gerrymander.Plot.Plotter(partitioner)

    # Animate
    Gerrymander.Plot.animate!(plotter)

end

main(ARGS)
