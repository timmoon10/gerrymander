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
const num_partitions::Int64 = 2
const relaxation_steps::Int64 = 100
const partisan_attraction::Float64 = 0.0
const partisan_repulsion::Float64 = 0.0
const personal_stdev::Float64 = 50.0 # km
const max_dist::Float64 = 200.0 # km
const bin_size::Float64 = 50.0 # km
const earth_radius::Float64 = 6370.286 # km (computed at (37.411764°N, 92.394544°W))
