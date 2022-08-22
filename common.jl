# Paths
project_dir = dirname(realpath(@__FILE__))
download_dir = joinpath(project_dir, "data", "downloads")
results_dir = joinpath(project_dir, "results")
county_data_file = joinpath(results_dir, "county_data.tsv")
geography_data_file = joinpath(results_dir, "geography.bin")
partition_file = joinpath(results_dir, "partition.tsv")
image_file = joinpath(results_dir, "partition.png")

### TODO Remove
results_dir            = project_dir * "/results/"
proto_file             = results_dir * "/gerrymander_pb.jl"
county_data_file       = results_dir * "/county_data.tsv"
interaction_graph_file = results_dir * "/interaction_graph.prototxt"
geography_graph_file   = results_dir * "/geography_graph.prototxt"
