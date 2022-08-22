#
# Construct geography graph.
# Nodes are geography regions and edges are border lengths.
#
import Serialization
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "graph.jl"))

# Import geography data
println("Importing geography data...")
regions = Serialization.deserialize(geography_data_file)

# Function to generate geography graph
function make_geography_graph(
    regions::Dict{Int64, Vector{Vector{Array{Float64, 2}}}}
    )::Dict{Int64, Dict{Int64, Float64}}

    # Construct empty graph
    graph = Dict{Int64, Dict{Int64, Float64}}(
        id => Dict{Int64, Float64}() for id in keys(regions))

    # Compute contribution of border segments to geography graph
    deg_to_rad::Float64 = pi / 180
    segment_owners = Dict{Tuple{Float64, Float64, Float64, Float64}, Int64}()
    for (id, region) in regions
        for polygon in region
            for border in polygon
                for i in 1:size(border,2)-1

                    # Border segment
                    x1 = border[1, i]
                    y1 = border[2, i]
                    x2 = border[1, i+1]
                    y2 = border[2, i+1]
                    segment = (x1, y1, x2, y2)

                    if haskey(segment_owners, segment)

                        # Another region shares this border segment
                        neighbor = segment_owners[segment]
                        if neighbor == id
                            continue
                        end
                        delete!(segment_owners, segment)

                        # Compute segment length with Haversine formula
                        l1 = deg_to_rad * x1
                        ph1 = deg_to_rad * y1
                        l2 = deg_to_rad * x2
                        ph2 = deg_to_rad * y2
                        hav_angle = (
                            sin((ph2-ph1)/2)^2
                            + cos(ph1) * cos(ph2) * sin((l2-l1)/2)^2
                        )
                        d = 2 * earth_radius * asin(sqrt(hav_angle))

                        # Add segment length to graph
                        if !haskey(graph[id], neighbor)
                            graph[id][neighbor] = 0.0
                            graph[neighbor][id] = 0.0
                        end
                        graph[id][neighbor] += d
                        graph[neighbor][id] = graph[id][neighbor]

                    else
                        # Wait for another region to encounter segment
                        segment_owners[(x2, y2, x1, y1)] = id
                    end
                end
            end
        end
    end

    return graph
end

# Construct geography graph
println("Constructing geography graph...")
graph = make_geography_graph(regions)

# Remove unconnected components of graph
println("Removing unconnected components of geography graph...")
search_start = minimum(keys(graph))
is_connected = find_connected(graph, search_start)
if !all(values(is_connected))
    println("Geography graph is unconnected!")
    graph = Dict{Int64, Dict{Int64, Float64}}(
        id => neighbors
        for (id, neighbors) in graph
        if is_connected[id]
    )
end

# Write results to file
println("Exporting geography graph...")
Serialization.serialize(geography_graph_file, graph)
