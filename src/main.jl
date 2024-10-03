include(joinpath(@__DIR__, "Gerrymander.jl"))

code = 10
println("State ", code, " is ", Gerrymander.DataFiles.state_names()[code])
