include(joinpath(@__DIR__, "Gerrymander.jl"))

function main()

    # State name
    state_id::UInt = 34
    state_name = Gerrymander.DataFiles.state_names()[state_id]
    println("State ", state_id, " is ", state_name)

    # Initial partition
    county_population_data = Gerrymander.DataFiles.load_county_populations([state_id])
    county_ids = Vector{UInt}(county_population_data[:, 1])
    partition = Dict(id => UInt(0) for id in county_ids)

    # Plot
    plotter = Gerrymander.Plot.Plotter(partition)
    Gerrymander.Plot.plot(plotter)

end

main()
