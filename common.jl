import ArgParse

# Command-line arguments
function parse_args()
    s = ArgParse.ArgParseSettings()
    @ArgParse.add_arg_table s begin
        "--num-partitions"
            arg_type = Int64
            default = 6
        "--relaxation-steps"
            arg_type = Int64
            default = 100
        "--partisan-attraction"
            arg_type = Float64
            default = 0.0
        "--partisan-repulsion"
            arg_type = Float64
            default = 0.0
        "--personal-stdev"
            arg_type = Float64
            default = 100.0 # km
    end
    return ArgParse.parse_args(s)
end
args = parse_args()

# Paths
project_dir = dirname(realpath(@__FILE__))
download_dir = joinpath(project_dir, "data", "downloads")
results_dir = joinpath(project_dir, "results")
county_data_file = joinpath(results_dir, "county_data.tsv")
geography_data_file = joinpath(results_dir, "geography.bin")
geography_graph_file = joinpath(results_dir, "geography_graph.bin")
interaction_graph_file = joinpath(results_dir, "interaction_graph.bin")
partition_file = joinpath(results_dir, "partition.tsv")
image_file = joinpath(results_dir, "partition.png")

# Constants
const num_partitions::Int64 = args["num-partitions"]
const relaxation_steps::Int64 = args["relaxation-steps"]
const partisan_attraction::Float64 = args["partisan-attraction"]
const partisan_repulsion::Float64 = args["partisan-repulsion"]
const personal_stdev::Float64 = args["personal-stdev"]
const max_dist::Float64 = 400.0 # km
const bin_size::Float64 = 100.0 # km
const earth_radius::Float64 = 6370.286 # km (computed at (37.411764°N, 92.394544°W))
