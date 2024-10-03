include(joinpath(@__DIR__, "Gerrymander.jl"))

state_id::UInt = 10
state_name = Gerrymander.DataFiles.state_names()[state_id]
println("State ", state_id, " is ", state_name)
Gerrymander.DataFiles.maybe_parse_county_populations([state_id])
