module Graph

import DataStructures
import Random

# Weighted graph
WeightedGraph = Dict{UInt, Dict{UInt, Float64}}

function make_partition_connected!(
    county_to_partition::Dict{UInt, UInt},
    partition_to_counties::Dict{UInt, Set{UInt}},
    county_populations::Dict{UInt, UInt},
    graph::WeightedGraph,
    )::Nothing

    # Iterate through partitions
    partition_ids = collect(keys(partition_to_counties))
    Random.shuffle!(partition_ids)
    for partition_id in partition_ids
        partition_counties = partition_to_counties[partition_id]

        # Split partition into connected regions with BFS
        connected_regions = Vector{Tuple{UInt, Set{UInt}}}()
        remaining_counties = copy(partition_counties)
        while !isempty(remaining_counties)
            start_county = pop!(remaining_counties)
            region = Set{UInt}([start_county])
            pop = county_populations[start_county]
            search_queue = DataStructures.Queue{UInt}()
            DataStructures.enqueue!(search_queue, start_county)
            while !isempty(search_queue)
                current = DataStructures.dequeue!(search_queue)
                for neighbor in keys(graph[current])
                    if in(neighbor, remaining_counties) && !in(neighbor, region)
                        push!(region, neighbor)
                        delete!(remaining_counties, neighbor)
                        pop += county_populations[neighbor]
                        DataStructures.enqueue!(search_queue, neighbor)
                    end
                end
            end
            push!(connected_regions, (pop, region))
        end

        # Nothing to be done if partition is already connected
        if length(connected_regions) <= 1
            continue
        end

        # Identify counties that disconnected from largest region
        (max_pop, max_region) = connected_regions[1]
        for (pop, region) in connected_regions[2:end]
            if pop > max_pop
                max_pop = pop
                max_region = region
            end
        end
        disconnected = setdiff(partition_counties, max_region)
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
                    delete!(disconnected, county_id)
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
