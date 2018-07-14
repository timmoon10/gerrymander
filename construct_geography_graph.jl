#!/usr/bin/julia
#
# Construct geography graph.
# Nodes are geography regions and edges are border lengths.
#
using ProtoBuf
import DataStructures
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Import geography data
println("Importing geography data...")
region_list_proto = MultiPolygonList()
open(geography_data_file, "r") do f
    readproto(f, region_list_proto)
end

# Function to compute border lengths
function update_border_lengths(id, border_proto,
                               segment_owners, border_lengths)
    const deg_to_rad = pi / 180
    num_coords = length(border_proto.x) - 1
    for i in 1:num_coords
        x1 = border_proto.x[i]
        y1 = border_proto.y[i]
        x2 = border_proto.x[i+1]
        y2 = border_proto.y[i+1]

        # Check if this segment has already been encountered
        # Note: If the segment has already been encountered, compute
        # the segment length. Otherwise record the segment and wait
        # for another region to encounter it.
        if haskey(segment_owners, (x1, y1, x2, y2))
            # Haversine formula for great circle distance
            neighbor = segment_owners[(x1, y1, x2, y2)]
            l1  = deg_to_rad * x1
            ph1 = deg_to_rad * y1
            l2  = deg_to_rad * x2
            ph2 = deg_to_rad * y2
            hav_angle = (sin((ph2-ph1)/2)^2
                         + cos(ph1) * cos(ph2) * sin((l2-l1)/2)^2)
            d = 2 * earth_radius * asin(sqrt(hav_angle))
            if haskey(border_lengths, neighbor)
                border_lengths[neighbor] += d
            else
                border_lengths[neighbor] = d
            end
        else
            segment_owners[(x2, y2, x1, y1)] = id
        end
    end
    return (segment_owners, border_lengths)
end    

# Construct geography graph
println("Constructing geography graph...")
graph = Dict{Int64, Dict{Int64, Float64}}()
id_list = Vector{Int64}()
segment_owners = Dict{Tuple{Float64, Float64, Float64, Float64}, Int64}()
for region_proto in region_list_proto.multi_polygon
    id = region_proto.id
    push!(id_list, id)

    # Compute border lengths
    border_lengths = Dict{Int64, Float64}()
    for polygon_proto in region_proto.polygon
        (segment_owners, border_lengths) = update_border_lengths(id,
                                                                 polygon_proto.exterior_border,
                                                                 segment_owners,
                                                                 border_lengths)
        for border_proto in polygon_proto.interior_border
            (segment_owners, border_lengths) = update_border_lengths(id,
                                                                     border_proto,
                                                                     segment_owners,
                                                                     border_lengths)
            
        end
    end

    # Construct graph edges corresponding to current region
    graph[id] = border_lengths
    for neighbor in keys(border_lengths)
        graph[neighbor][id] = border_lengths[neighbor]
    end
    
end

# BFS to find connected component of graph
println("Breadth-first search on geography graph...")
search_start = id_list[end]
search_queue = DataStructures.list(search_start)
is_connected = Dict{Int64, Bool}()
for id in id_list
    is_connected[id] = false
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
println("Converting geography graph to protobuf...")
graph_proto = Graph()
fillset(graph_proto, :node)
fillset(graph_proto, :edge)
for id in id_list
    if is_connected[id]
        add_field!(graph_proto, :node, id)
        weights = graph[id]
        for neighbor in keys(weights)
            if is_connected[neighbor] && id < neighbor
                edge_proto = GraphEdge()
                set_field!(edge_proto, :node1, id)
                set_field!(edge_proto, :node2, neighbor)
                set_field!(edge_proto, :weight, weights[neighbor])
                add_field!(graph_proto, :edge, edge_proto)
            end
        end
    end
end

# Write results to file
println("Exporting geography graph...")
open(geography_graph_file, "w") do f
    writeproto(f, graph_proto)
end
