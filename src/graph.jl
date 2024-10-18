module Graph

import DataStructures
import Random

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

function random_partition(
    num_partitions::UInt,
    graph::Graph.WeightedGraph,
    )::Tuple{Dict{UInt, UInt}, Dict{UInt, Set{UInt}}}

    # Check number of partitions
    county_ids = collect(keys(graph))
    num_counties = length(county_ids)
    if num_partitions > num_counties || num_partitions == 0
        throw(
            "Invalid number of partitions "
            * "($num_partitions partitions, $num_counties counties)"
        )
    end

    # Partition data
    partition_ids = collect(1:num_partitions)
    county_to_partition = Dict{UInt, UInt}()
    partition_to_counties = Dict{UInt, Set{UInt}}(
        id => Set{UInt}() for id in partition_ids)

    # Objects for tracking partition neighbors
    assigned_counties = Set{UInt}()
    sizehint!(assigned_counties, num_counties)
    partition_edges = Dict{UInt, Set{UInt}}(
        id => Set{UInt}() for id in partition_ids)
    unfinished_partitions = Set{UInt}()

    # Seed partitions with randomly chosen counties
    seed_counties = [
        county_ids[i] for i in Random.randperm(num_counties)[1:num_partitions]]
    union!(assigned_counties, seed_counties)
    for (partition_id, county_id) in enumerate(seed_counties)
        county_to_partition[county_id] = partition_id
        push!(partition_to_counties[partition_id], county_id)
        for neighbor_id in keys(graph[county_id])
            if !in(neighbor_id, assigned_counties)
                push!(partition_edges[partition_id], county_id)
                push!(unfinished_partitions, partition_id)
                break
            end
        end
    end

    function unassigned_neighbors(county_id::UInt)::Vector{UInt}
        out = Vector{UInt}()
        for neighbor_id in keys(graph[county_id])
            if !in(neighbor_id, assigned_counties)
                push!(out, neighbor_id)
            end
        end
        return out
    end

    function has_unassigned_neighbors(county_id::UInt)::Bool
        for neighbor_id in keys(graph[county_id])
            if !in(neighbor_id, assigned_counties)
                return true
            end
        end
        return false
    end

    while !isempty(unfinished_partitions)

        # Add random county to random partition
        partition_id = rand(unfinished_partitions)
        edge_county_id = rand(partition_edges[partition_id])
        county_id = rand(unassigned_neighbors(edge_county_id))
        county_to_partition[county_id] = partition_id
        push!(partition_to_counties[partition_id], county_id)

        # Check if county is at partition edge
        push!(assigned_counties, county_id)
        if has_unassigned_neighbors(county_id)
            push!(partition_edges[partition_id], county_id)
        end

        # Update whether neighbors are still at partition edge
        for neighbor_id in keys(graph[county_id])
            if !in(neighbor_id, assigned_counties)
                continue
            end
            neighbor_partition_id = county_to_partition[neighbor_id]
            if !in(neighbor_id, partition_edges[neighbor_partition_id])
                continue
            end
            if has_unassigned_neighbors(neighbor_id)
                continue
            end
            pop!(partition_edges[neighbor_partition_id], neighbor_id)
            if isempty(partition_edges[neighbor_partition_id])
                pop!(unfinished_partitions, neighbor_partition_id)
            end
        end

    end

    # Check results
    if length(county_to_partition) != num_counties
        throw("Failed to partition graph, maybe because graph is disconnected")
    end

    return (county_to_partition, partition_to_counties)

end

end  # module Graph
