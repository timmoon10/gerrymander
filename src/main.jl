include(joinpath(@__DIR__, "Gerrymander.jl"))

function main()

    # State name
    state_id::UInt = 34
    state_name = Gerrymander.DataFiles.state_names()[state_id]
    println("State ", state_id, " is ", state_name)
    state_ids = [state_id]

    # Partitioner
    num_partitions::UInt = 3
    personal_stdev::Float64 = 100.0
    max_distance::Float64 = 300.0
    partitioner = Gerrymander.SimulatedAnnealing.Partitioner(
        num_partitions,
        state_ids,
        personal_stdev,
        max_distance,
    )

    # Plotter
    plotter = Gerrymander.Plot.Plotter(partitioner)

    # Perform simulated annealing
    num_steps = 100
    for step in 1:num_steps
        update_plotter = step == num_steps
        Gerrymander.SimulatedAnnealing.step!(partitioner)
    end

    # Plot
    Gerrymander.Plot.plot(plotter)

end

main()
