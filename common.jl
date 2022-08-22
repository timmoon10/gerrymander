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
const earth_radius::Float64 = 6370.286 # km (computed at (37.411764°N, 92.394544°W))
