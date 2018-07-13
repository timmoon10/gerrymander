#!/usr/bin/julia

# Project parameters
const num_partitions      = 6
const balance_tolerance   = 2.0
const partisan_attraction = 4.0
const partisan_repulsion  = 4.0
const personal_stdev      = 50.0     # km
const max_dist            = 200.0    # km
const bin_size            = 25.0     # km
const earth_radius        = 6370.286 # km (computed at (37.411764°N, 92.394544°W))

# Project files
project_dir            = dirname(@__FILE__) * "/"
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
