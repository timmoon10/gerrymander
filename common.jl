#!/usr/bin/julia

# Project parameters
const num_partitions      = 6
const balance_tolerance   = 1.5
const partisan_attraction = 4.0
const partisan_repulsion  = 8.0
const personal_stdev      = 25.0     # km
const max_dist            = 200.0    # km
const bin_size            = 25.0     # km
const earth_radius        = 6370.286 # km (computed at (37.411764°N, 92.394544°W))

# Project files
project_dir            = AbstractString(dirname(@__FILE__)) * "/"
download_dir           = project_dir * "/data/downloads/"
results_dir            = project_dir * "/results/"
proto_file             = results_dir * "/gerrymander_pb.jl"
county_data_file       = results_dir * "/county_data.tsv"
county_bounds_file     = results_dir * "/bounds.prototxt"
interaction_graph_file = results_dir * "/interaction_graph.prototxt"
geography_graph_file   = results_dir * "/geography_graph.prototxt"
partition_file         = results_dir * "/partition.tsv"
image_file             = results_dir * "/parts.png"

# External dependencies
protoc                 = "/home/moon/src/protobuf/src/protoc"
julia_protobuf_dir     = "/home/moon/.julia/v0.4/ProtoBuf"
color_generator        = "/home/moon/src/randomcolor-py/getcolor.py"

# Create directories if needed
if !isdir(download_dir)
    mkdir(download_dir)
end
if !isdir(results_dir)
    mkdir(results_dir)
end

# Initialize protobuf
println("Initializing protobuf...")
if !isfile(proto_file)
    run(`$protoc --plugin=$julia_protobuf_dir/plugin/protoc-gen-julia 
         -I=$project_dir --julia_out=$proto_file
         $project_dir/gerrymander.proto`)
end
include(proto_file)
