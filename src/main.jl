include(joinpath(@__DIR__, "Gerrymander.jl"))

function main()
    state_id::UInt = 34
    state_name = Gerrymander.DataFiles.state_names()[state_id]
    println("State ", state_id, " is ", state_name)
    Gerrymander.DataFiles.load_county_populations([state_id])
    Gerrymander.DataFiles.load_county_boundaries()
end

main()
