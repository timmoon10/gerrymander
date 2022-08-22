import DataStructures

# Function to find connected component of graph with BFS
function find_connected(
    graph::Dict{Int64, Dict{Int64, Float64}},
    search_start::Int64,
    )::Dict{Int64, Bool}
    is_connected = Dict{Int64, Bool}(id => false for id in keys(graph))
    is_connected[search_start] = true
    search_queue = DataStructures.Queue{Int64}()
    DataStructures.enqueue!(search_queue, search_start)
    while !isempty(search_queue)
        current = DataStructures.dequeue!(search_queue)
        for neighbor in keys(graph[current])
            if !is_connected[neighbor]
                DataStructures.enqueue!(search_queue, neighbor)
                is_connected[neighbor] = true
            end
        end
    end
    return is_connected
end
