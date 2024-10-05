include(joinpath(@__DIR__, "Gerrymander.jl"))

function main()
    state_id::UInt = 34
    state_name = Gerrymander.DataFiles.state_names()[state_id]
    println("State ", state_id, " is ", state_name)
    county_population_data = Gerrymander.DataFiles.load_county_populations([state_id])
    county_ids = Vector{UInt}(county_population_data[:, 1])
    plotter = Gerrymander.Plot.Plotter(county_ids)
    Gerrymander.Plot.plot(plotter)
end

main()
