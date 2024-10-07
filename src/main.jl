include(joinpath(@__DIR__, "Gerrymander.jl"))

function main()

    # State name
    state_id::UInt = 34
    state_name = Gerrymander.DataFiles.state_names()[state_id]
    println("State ", state_id, " is ", state_name)
    state_ids = [state_id]

    # Graphs
    county_adjacency_graph = Gerrymander.DataGraphs.county_adjacency_graph(state_ids)
    county_interaction_graph = Gerrymander.DataGraphs.county_interaction_graph(
        state_ids,
        100.0,  # personal_stdev
        300.0,  # max_distance
    )

    # Initial partition
    county_population_data = Gerrymander.DataFiles.load_county_populations(state_ids)
    county_ids = Vector{UInt}(county_population_data[:, 1])
    partition = Dict(id => UInt(0) for id in county_ids)

    # Plot
    plotter = Gerrymander.Plot.Plotter(partition)
    Gerrymander.Plot.plot(plotter)

end

main()
