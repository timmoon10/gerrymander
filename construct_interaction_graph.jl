#!/usr/bin/julia
#
# Construct county interaction matrix
#
using Iterators, ProtoBuf
import DataStructures
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Import county population data
println("Importing county population data...")
(county_data, _) = readdlm(county_data_file, '\t', header=true)

# Construct graph one county at a time
println("Constructing county interaction graph...")
graph = Dict{Int64, Dict{Int64, Float64}}()
geoid_list = DataStructures.list()
binned_county_data = Dict{Tuple{Int64, Int64, Int64},
                          Vector{Tuple{Int64, Int64, Int64, Int64,
                                       Float64, Float64, Float64, Float64}}}()
const personal_var = personal_stdev^2
const max_dist2 = max_dist^2
for row in 1:size(county_data, 1)
    edge_weights = Dict{Int64, Float64}()

    # Current county data
    geoid1     = Int64(county_data[row, 1])
    pop1       = Int64(county_data[row, 2])
    dem_votes1 = Int64(county_data[row, 3])
    gop_votes1 = Int64(county_data[row, 4])
    x1         = Float64(county_data[row, 5])
    y1         = Float64(county_data[row, 6])
    z1         = Float64(county_data[row, 7])
    var1       = Float64(county_data[row, 8])

    # Compute edge weights for counties in nearby bins
    pos1 = [x1, y1, z1]
    bins_min = floor(Int64, (pos1 - max_dist) / bin_size)
    bins_max = floor(Int64, (pos1 + max_dist) / bin_size)
    for bin in Iterators.product(bins_min[1]:bins_max[1],
                                 bins_min[2]:bins_max[2],
                                 bins_min[3]:bins_max[3])
        if !haskey(binned_county_data, bin)
            continue
        end
        for (geoid2, pop2, dem_votes2, gop_votes2, x2, y2, z2, var2) in binned_county_data[bin]
            dist2 = (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2
            if dist2 > max_dist2
                continue
            end
            interaction_var = var1 + var2 + 2 * personal_var
            unit_interaction = (exp(-0.5 * dist2 / interaction_var)
                                / (2 * pi * sqrt(interaction_var)))
            eff_pop_product = (pop1 * pop2
                               + partisan_attraction * (dem_votes1 * dem_votes2
                                                        + gop_votes1 * gop_votes2)
                               - partisan_repulsion * (dem_votes1 * gop_votes2
                                                       + gop_votes1 * dem_votes2))
            edge_weights[geoid2] = eff_pop_product * unit_interaction
        end
    end

    # Add current county data to graph and binned cache
    graph[geoid1] = edge_weights
    for geoid2 in keys(edge_weights)
        graph[geoid2][geoid1] = edge_weights[geoid2]
    end
    geoid_list = DataStructures.cons(geoid1, geoid_list)
    bin = (floor(Int64, x1 / bin_size),
           floor(Int64, y1 / bin_size),
           floor(Int64, z1 / bin_size))
    if !haskey(binned_county_data, bin)
        binned_county_data[bin] = Vector{Tuple{Int64, Int64, Int64, Int64,
                                               Float64, Float64, Float64, Float64}}(0)
    end
    push!(binned_county_data[bin],
          (geoid1, pop1, dem_votes1, gop_votes1, x1, y1, z1, var1))
    
end

# BFS to find connected component of graph
println("Breadth-first search on county interaction graph...")
search_start = geoid_list.head
search_queue = DataStructures.list(search_start)
is_connected = Dict{Int64, Bool}()
for geoid in geoid_list
    is_connected[geoid] = false
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
println("Converting county interaction graph to protobuf...")
graph_proto = Graph()
fillset(graph_proto, :node)
fillset(graph_proto, :edge)
for geoid in geoid_list
    if is_connected[geoid]
        add_field!(graph_proto, :node, geoid)
        neighbor_list = graph[geoid]
        for neighbor in keys(neighbor_list)
            if is_connected[neighbor] && geoid < neighbor
                edge_proto = GraphEdge()
                set_field!(edge_proto, :node1, geoid)
                set_field!(edge_proto, :node2, neighbor)
                set_field!(edge_proto, :weight, neighbor_list[neighbor])
                add_field!(graph_proto, :edge, edge_proto)
            end
        end
    end
end

# Write results to file
println("Writing county interaction graph to file...")
open(interaction_graph_file, "w") do f
    writeproto(f, graph_proto)
end
