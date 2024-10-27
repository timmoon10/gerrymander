include(joinpath(dirname(@__DIR__), "src", "Gerrymander.jl"))
import ArgParse

function parse_args()
    settings = ArgParse.ArgParseSettings()
    @ArgParse.add_arg_table settings begin
        "--num-partitions"
            arg_type = UInt
            default = 4
        "--personal-stdev"
            arg_type = Float64
            default = 100.0
        "--personal-max-distance"
            arg_type = Float64
            default = 400.0
        "states"
            nargs = '*'
            default = ["6"]  # California
    end
    return ArgParse.parse_args(settings)
end

function main()

    # Parse command-line arguments
    args = parse_args()

    # Parse states
    state_ids = Set{UInt}()
    state_names = Gerrymander.DataFiles.state_names()
    if length(args["states"]) == 1 && args["states"][1] == "all"
        union!(state_ids, keys(state_names))
    else
        for arg in args["states"]
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
    println()
    println(message)

    # Partitioner
    partitioner = Gerrymander.SimulatedAnnealing.Partitioner(
        args["num-partitions"],
        state_ids,
        args["personal-stdev"],
        args["personal-max-distance"],
    )

    # Plotter
    plotter = Gerrymander.Plot.Plotter(partitioner)

    # Animate
    Gerrymander.Plot.animate!(plotter)

end

# Run main function
main()
