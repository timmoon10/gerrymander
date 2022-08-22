import DataStructures

# Function to find connected component of graph with BFS
function find_connected_vertices(
    graph::Dict{Int64, Dict{Int64, Float64}},
    search_start::Int64,
    )::Set{Int64}
    found_set = Set{Int64}()
    search_queue = DataStructures.Queue{Int64}()
    push!(found_set, search_start)
    DataStructures.enqueue!(search_queue, search_start)
    while !isempty(search_queue)
        current = DataStructures.dequeue!(search_queue)
        for neighbor in keys(graph[current])
            if !in(neighbor, found_set)
                DataStructures.enqueue!(search_queue, neighbor)
                push!(found_set, neighbor)
            end
        end
    end
    return found_set
end

function construct_subgraph(
    graph::Dict{Int64, Dict{Int64, Float64}},
    vertices::Set{Int64},
    )::Dict{Int64, Dict{Int64, Float64}}
    return Dict{Int64, Dict{Int64, Float64}}(
        id => Dict{Int64, Float64}(
            neighbor => weight
            for (neighbor, weight) in graph[id]
            if in(neighbor, vertices)
        )
        for id in vertices
    )
end
