module Graph

import DataStructures

# Weighted graph
WeightedGraph = Dict{UInt, Dict{UInt, Float64}}

function make_partition_contiguous!(
    county_to_partition::Dict{UInt, UInt},
    partition_to_counties::Dict{UInt, Set{UInt}},
    partition_center_county::Dict{UInt, UInt},
    graph::WeightedGraph,
    )::Nothing

    # Iterate through partitions
    for partition_id in keys(partition_to_counties)
        partition_counties = partition_to_counties[partition_id]

        # BFS from partition center
        center_id = partition_center_county[partition_id]
        found_set = Set{UInt}([center_id])
        search_queue = DataStructures.Queue{UInt}()
        DataStructures.enqueue!(search_queue, center_id)
        while !isempty(search_queue)
            current = DataStructures.dequeue!(search_queue)
            for neighbor in keys(graph[current])
                if in(neighbor, partition_counties) && !in(neighbor, found_set)
                    push!(found_set, neighbor)
                    DataStructures.enqueue!(search_queue, neighbor)
                end
            end
        end

        # Identify disconnected counties
        disconnected = setdiff(partition_counties, found_set)
        setdiff!(partition_counties, disconnected)

        # Reassign disconnected counties to neighboring partitions
        while !isempty(disconnected)
            for county_id in collect(disconnected)
                for neighbor in keys(graph[county_id])
                    neighbor_partition_id = county_to_partition[neighbor]
                    if neighbor_partition_id == partition_id
                        continue
                    end
                    county_to_partition[county_id] = neighbor_partition_id
                    push!(
                        partition_to_counties[neighbor_partition_id],
                        county_id,
                    )
                    break
                end
            end
        end

    end

end

end  # module Graph
