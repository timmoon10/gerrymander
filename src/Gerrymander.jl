module Gerrymander

include(joinpath(@__DIR__, "constants.jl"))
include(joinpath(@__DIR__, "graph.jl"))
include(joinpath(@__DIR__, "geometry.jl"))
include(joinpath(@__DIR__, "data_files.jl"))
include(joinpath(@__DIR__, "data_graphs.jl"))
include(joinpath(@__DIR__, "simulated_annealing.jl"))
using .SimulatedAnnealing: step!
include(joinpath(@__DIR__, "plot.jl"))

end  # module Gerrymander
