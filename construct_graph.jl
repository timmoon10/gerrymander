#!/usr/bin/julia

#
# Construct adjacency matrix
#
# The radius of the Earth is computed with
# https://rechneronline.de/earth-radius/ at the population center of
# the U.S. (37.411764°N, 92.394544°W).

# Parameters
project_dir  = "/home/moon/Documents/gerrymander/"
output_dir   = project_dir * "/output"
bounds_file  = output_dir * "/bounds.prototxt"
graph_file   = output_dir * "/graph.prototxt"
earth_radius = 6370286 # meters

# Import packages
using ProtoBuf
import DataStructures

# Read boundary data from protobuf
println("Reading boundary data...")
include(output_dir * "/gerrymander_pb.jl")
bounds_list_proto = CountyBoundariesList()
open(bounds_file, "r") do f
    readproto(f, bounds_list_proto)
end

# Construct county graph one county at a time
println("Constructing county graph...")
graph = Dict{UInt32, Dict{UInt32, Float32}}()
county_list = DataStructures.list()
segment_list = Dict{Tuple{Float32, Float32, Float32, Float32}, UInt32}()
deg_to_rad = pi / 180
for bounds_proto in bounds_list_proto.county_bounds
    county = bounds_proto.geoid
    county_list = DataStructures.cons(county, county_list)
    neighbor_list = Dict{UInt32, Float32}()

    for polygon_proto in bounds_proto.polygon
        num_coords = length(polygon_proto.x) - 1
        for i = 1:num_coords

            # Check if another county has encountered this segment
            # Note: If another county has seen it, then compute the
            # segment length. Otherwise record the segment and wait
            # for another county to encounter it.
            x1 = polygon_proto.x[i]
            y1 = polygon_proto.y[i]
            x2 = polygon_proto.x[i+1]
            y2 = polygon_proto.y[i+1]
            if haskey(segment_list, (x1, y1, x2, y2))
                # Compute length of boundary segment and record
                # Note: haversine formula for great circle distance
                neighbor = segment_list[(x1, y1, x2, y2)]
                l1  = deg_to_rad * x1
                ph1 = deg_to_rad * y1
                l2  = deg_to_rad * x2
                ph2 = deg_to_rad * y2
                hav_angle = (sin((ph2-ph1)/2)^2
                             + cos(ph1) * cos(ph2) * sin((l2-l1)/2)^2)
                d = 2 * earth_radius * asin(sqrt(hav_angle))
                if haskey(neighbor_list, neighbor)
                    neighbor_list[neighbor] += d
                else
                    neighbor_list[neighbor] = d
                end
            else
                # Add segment to list of encountered segments
                segment_list[(x2, y2, x1, y1)] = county
            end
        end
    end

    # Construct graph edges corresponding to current county
    graph[county] = neighbor_list
    for neighbor in keys(neighbor_list)
        graph[neighbor][county] = neighbor_list[neighbor]
    end
    
end

# BFS to find connected component of graph
println("Breadth-first search on county graph...")
search_start = county_list.head
search_queue = DataStructures.list(search_start)
is_connected = Dict{UInt32, Bool}()
for county in county_list
    is_connected[county] = false
end
is_connected[search_start] = true
while !isempty(search_queue)
    current = search_queue.head
    search_queue = search_queue.tail
    for neighbor in keys(graph[current])
        if !is_connected[neighbor]
            search_queue = DataStructures.cons(neighbor, search_queue)
            is_connected[neighbor] = true
        end
    end
end

# Convert graph to protobuf format
println("Converting county graph to protobuf...")
graph_proto = Graph()
fillset(graph_proto, :node)
fillset(graph_proto, :edge)
for county in county_list
    if is_connected[county]
        add_field!(graph_proto, :node, county)
        neighbor_list = graph[county]
        for neighbor in keys(neighbor_list)
            if is_connected[neighbor] && county < neighbor
                edge_proto = GraphEdge()
                set_field!(edge_proto, :node1, county)
                set_field!(edge_proto, :node2, neighbor)
                set_field!(edge_proto, :weight, neighbor_list[neighbor])
                add_field!(graph_proto, :edge, edge_proto)
            end
        end
    end
end

# Write results to file
println("Writing county graph to file...")
open(graph_file, "w") do f
    writeproto(f, graph_proto)
end
