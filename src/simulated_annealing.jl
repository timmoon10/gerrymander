module SimulatedAnnealing

import Base.Threads
import DataStructures
import Random

import ..DataFiles
import ..DataGraphs
import ..Graph

function contiguous_partition!(
    county_to_partition::Dict{UInt, UInt},
    partition_to_counties::Dict{UInt, Set{UInt}},
    partition_center_county::Dict{UInt, UInt},
    graph::Graph.WeightedGraph,
    )

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

function voronoi_partition(
    num_partitions::UInt,
    graph::Graph.WeightedGraph,
    state_ids::AbstractVector{UInt},
    )::Tuple{Dict{UInt, UInt}, Dict{UInt, Set{UInt}}}

    # County data
    county_population_data = DataFiles.load_county_populations(state_ids)
    num_counties = size(county_population_data, 1)
    county_ids = Vector{UInt}(county_population_data[:, 1])
    county_points = Dict{UInt, Vector{Float64}}(
        county_ids[i] => county_population_data[i, 3:5]
        for i in 1:num_counties)

    # Partition data
    partition_ids = collect(1:num_partitions)
    county_to_partition = Dict{UInt, UInt}()
    partition_to_counties = Dict{UInt, Set{UInt}}(
        id => Set{UInt}() for id in partition_ids)

    "Square of Euclidean distance between 3D points"
    function dist_sq(
        point1::Vector{Float64},
        point2::Vector{Float64},
        )::Float64
        @inbounds dx = point1[1] - point2[1]
        @inbounds dy = point1[2] - point2[2]
        @inbounds dz = point1[3] - point2[3]
        return dx * dx + dy * dy + dz * dz
    end

    # Randomly pick counties
    center_ids = [
        county_ids[i]
        for i in Random.randperm(num_counties)[1:num_partitions]]
    center_points = [county_points[id] for id in center_ids]
    for (partition_id, county_id) in enumerate(center_ids)
        county_to_partition[county_id] = partition_id
        push!(partition_to_counties[partition_id], county_id)
    end

    # Assign points by Voronoi region
    for (county_id, point) in county_points
        if haskey(county_to_partition, county_id)
            continue
        end
        partition_id = argmin(
            [dist_sq(point, center) for center in center_points])
        county_to_partition[county_id] = partition_id
        push!(partition_to_counties[partition_id], county_id)
    end

    # Make sure partitions are contiguous
    graph = DataGraphs.county_adjacency_graph(state_ids)
    center_ids = Dict{UInt, UInt}(
        partition_id => county_id
        for (partition_id, county_id) in enumerate(center_ids))
    contiguous_partition!(
        county_to_partition,
        partition_to_counties,
        center_ids,
        graph,
    )

    return (county_to_partition, partition_to_counties)

end

mutable struct Partitioner
    temperature::Float64
    population_weight::Float64
    adjacency_graph::Graph.WeightedGraph
    interaction_graph::Graph.WeightedGraph
    county_to_partition::Dict{UInt, UInt}
    partition_to_counties::Dict{UInt, Set{UInt}}
    county_populations::Dict{UInt, UInt}
    partition_populations::Dict{UInt, UInt}
    swap_candidates::Dict{UInt, Dict{UInt, Float64}}
    on_key_func::Union{Function, Nothing}
end

function Partitioner(
    num_partitions::UInt,
    state_ids::AbstractVector{UInt},
    interaction_personal_stdev::Float64,
    interaction_max_distance::Float64;
    temperature::Float64 = 1.0,
    population_weight::Float64 = 1.0,
    )::Partitioner

    # Data graphs
    adjacency_graph = DataGraphs.county_adjacency_graph(state_ids)
    interaction_graph = DataGraphs.county_interaction_graph(
        state_ids,
        interaction_personal_stdev,
        interaction_max_distance,
    )

    # County data
    county_population_data = DataFiles.load_county_populations(state_ids)
    num_counties = size(county_population_data, 1)
    county_ids = Vector{UInt}(county_population_data[:, 1])
    county_populations = Dict{UInt, UInt}(
        county_ids[i] => county_population_data[i, 2]
        for i in 1:num_counties)

    # Initial partition
    (county_to_partition, partition_to_counties) = voronoi_partition(
        num_partitions,
        adjacency_graph,
        state_ids,
    )
    partition_populations = Dict{UInt, UInt}()
    sizehint!(partition_populations, Int(num_partitions))
    for (partition_id, counties) in partition_to_counties
        pop::UInt = 0
        for county_id in counties
            pop += county_populations[county_id]
        end
        partition_populations[partition_id] = pop
    end

    # Candidate county swaps
    swaps_candidates = Dict{UInt, Dict{UInt, Float64}}(
        id => Dict{UInt, Float64}() for id in county_ids)

    # Construct partitioner object
    partitioner = Partitioner(
        temperature,
        population_weight,
        adjacency_graph,
        interaction_graph,
        county_to_partition,
        partition_to_counties,
        county_populations,
        partition_populations,
        swaps_candidates,
        nothing,
    )

    "Logic for key presses"
    function on_key(event)
        if event.key == "h"
            println("Commands")
            println("--------")
            println("h: help message")
            println("esc: exit")
            println("p: pause animation")
            println("i: partitioner state")
            println("-/=: adjust temperature")
            println("[/]: adjust population weight")
        elseif event.key == "i"
            print_info(partitioner)
        elseif event.key == "-"
            partitioner.temperature /= 2
            println("Temperature: ", partitioner.temperature)
        elseif event.key == "="
            partitioner.temperature *= 2
            println("Temperature: ", partitioner.temperature)
        elseif event.key == "["
            partitioner.population_weight /= 2
            println("Population weight: ", partitioner.population_weight)
        elseif event.key == "]"
            partitioner.population_weight *= 2
            println("Population weight: ", partitioner.population_weight)
        end
    end

    # Register logic for key presses
    partitioner.on_key_func = on_key

    # Find swap candidates
    @Base.Threads.threads for county_id in county_ids
        update_county_swap_candidates!(partitioner, county_id)
    end

    return partitioner

end

function print_info(partitioner::Partitioner)

    # Print partitioner state
    println("Partitioner properties")
    println("----------------------")
    println("Temperature: ", partitioner.temperature)
    println("Population weight: ", partitioner.population_weight)
    print("\n")

    # Print partition populations
    println("Partition populations")
    println("---------------------")
    partition_ids = collect(keys(partitioner.partition_populations))
    sort!(partition_ids)
    for partition_id in partition_ids
        pop::Int = partitioner.partition_populations[partition_id]
        println("Partition ", Int(partition_id), ": ", pop)
    end
    print("\n")

    # Print partition populations
    println("Partition affinites")
    println("---------------------")
    for partition_id in partition_ids
        partition_affinity::Float64 = 0
        partition_counties = partitioner.partition_to_counties[partition_id]
        for county_id in partition_counties
            for (neighbor_id, affinity) in partitioner.interaction_graph[county_id]
                if neighbor_id != county_id && in(neighbor_id, partition_counties)
                    partition_affinity += affinity
                end

            end
        end
        println("Partition ", Int(partition_id), ": ", partition_affinity)
    end
    print("\n")

end

function update_county_swap_candidates!(
    partitioner::Partitioner,
    county_id::UInt,
    )

    # Reset list of swap candidates
    swap_candidates = partitioner.swap_candidates[county_id]
    empty!(swap_candidates)

    # Get partitions adjacent to county
    partition_id = partitioner.county_to_partition[county_id]
    neighbor_partitions = Set{UInt}([partition_id])
    for (neighbor_id, _) in partitioner.adjacency_graph[county_id]
        push!(
            neighbor_partitions,
            partitioner.county_to_partition[neighbor_id],
        )
    end

    # Return immediately if county is within partition interior
    if length(neighbor_partitions) == 1
        return
    end

    # Compute affinity to partitions
    partition_affinities = DataStructures.OrderedDict{UInt, Float64}(
        id => 0 for id in neighbor_partitions)
    for (neighbor_id, affinity) in partitioner.interaction_graph[county_id]
        if neighbor_id == county_id
            continue
        end
        neighbor_partition_id = partitioner.county_to_partition[neighbor_id]
        if haskey(partition_affinities, neighbor_partition_id)
            partition_affinities[neighbor_partition_id] += affinity
        end
    end

    # Update list of swap candidates
    self_affinity = partition_affinities[partition_id]
    delete!(neighbor_partitions, partition_id)
    for id in neighbor_partitions
        swap_candidates[id] = partition_affinities[id] - self_affinity
    end

end

function can_swap_county(
    partitioner::Partitioner,
    county_id::UInt,
    partition_id::UInt,
    )::Bool
    src_partition_id = partitioner.county_to_partition[county_id]
    dst_partition_id = partition_id

    # Get neighboring counties
    neighbors = keys(partitioner.adjacency_graph[county_id])
    src_neighbors = Set{UInt}()
    dst_neighbors = Set{UInt}()
    sizehint!(src_neighbors, length(neighbors))
    sizehint!(dst_neighbors, length(neighbors))
    for neighbor_id in neighbors
        neighbor_partition_id = partitioner.county_to_partition[neighbor_id]
        if neighbor_partition_id == src_partition_id
            push!(src_neighbors, neighbor_id)
        elseif neighbor_partition_id == dst_partition_id
            push!(dst_neighbors, neighbor_id)
        end
    end

    # Cannot swap isolated counties
    if isempty(src_neighbors) || isempty(dst_neighbors)
        return false
    end

    # Local BFS to check if county is local cut vertex
    not_visited = src_neighbors
    start_id = pop!(not_visited)
    search_queue = DataStructures.Queue{UInt}()
    DataStructures.enqueue!(search_queue, start_id)
    while !isempty(search_queue)
        current = DataStructures.dequeue!(search_queue)
        for neighbor in keys(partitioner.adjacency_graph[current])
            if in(neighbor, not_visited)
                delete!(not_visited, neighbor)
                DataStructures.enqueue!(search_queue, neighbor)
            end
        end
    end
    if !isempty(not_visited)
        return false
    end

    # Local BFS to check if county would become local cut vertex
    not_visited = dst_neighbors
    start_id = pop!(not_visited)
    search_queue = DataStructures.Queue{UInt}()
    DataStructures.enqueue!(search_queue, start_id)
    while !isempty(search_queue)
        current = DataStructures.dequeue!(search_queue)
        for neighbor in keys(partitioner.adjacency_graph[current])
            if in(neighbor, not_visited)
                delete!(not_visited, neighbor)
                DataStructures.enqueue!(search_queue, neighbor)
            end
        end
    end
    if !isempty(not_visited)
        return false
    end

    # County swap is valid
    return true

end

function swap_county!(
    partitioner::Partitioner,
    county_id::UInt,
    partition_id::UInt;
    )

    # Source and destination partitions
    src_partition_id = partitioner.county_to_partition[county_id]
    dst_partition_id = partition_id
    if dst_partition_id == src_partition_id
        return
    end

    # Update partitions
    partitioner.county_to_partition[county_id] = dst_partition_id
    delete!(partitioner.partition_to_counties[src_partition_id], county_id)
    push!(partitioner.partition_to_counties[dst_partition_id], county_id)

    # Update partition populations
    county_population = partitioner.county_populations[county_id]
    partitioner.partition_populations[src_partition_id] -= county_population
    partitioner.partition_populations[dst_partition_id] += county_population

    # Update swap candidates
    neighborhood = Set{UInt}([county_id])
    for neighbor in keys(partitioner.adjacency_graph[county_id])
        push!(neighborhood, neighbor)
    end
    for neighbor in keys(partitioner.interaction_graph[county_id])
        push!(neighborhood, neighbor)
    end
    @Base.Threads.threads for county_id in collect(neighborhood)
        update_county_swap_candidates!(partitioner, county_id)
    end

end

function step!(partitioner::Partitioner)

    # Return immediately if there are no candidate swaps
    num_swap_candidates::Int = 0
    for county_swaps in values(partitioner.swap_candidates)
        num_swap_candidates += length(county_swaps)
    end
    if num_swap_candidates == 0
        return
    end

    # Compute population statistics
    total_population::UInt = 0
    sum_population_sq::Float64 = 0
    for pop in values(partitioner.partition_populations)
        total_population += pop
        pop_float::Float64 = pop
        sum_population_sq += pop_float * pop_float
    end
    num_partitions = length(partitioner.partition_populations)
    mean_population = Float64(total_population) / num_partitions
    mean_population_sq = sum_population_sq / num_partitions
    var_population = mean_population_sq - mean_population * mean_population
    var_population = max(var_population, 1)
    stdev_population = sqrt(var_population)

    # Population scores
    population_scores = Dict{UInt, Float64}()
    for (partition_id, pop) in partitioner.partition_populations
        pop_zscore = (pop - mean_population) / stdev_population
        population_scores[partition_id] = -partitioner.population_weight * pop_zscore
    end

    # Flatten candidate swaps and compute affinity statistics
    swaps = Vector{Tuple{UInt, UInt}}()
    scores = Vector{Float64}()
    sizehint!(swaps, num_swap_candidates)
    sizehint!(scores, num_swap_candidates)
    sum_affinity::Float64 = 0
    sum_affinity_sq::Float64 = 0
    for (county_id, county_swaps) in partitioner.swap_candidates
        for (partition_id, affinity) in county_swaps
            push!(swaps, (county_id, partition_id))
            push!(scores, affinity)
            sum_affinity += affinity
            sum_affinity_sq += affinity * affinity
        end
    end
    mean_affinity = sum_affinity / num_swap_candidates
    mean_affinity_sq = sum_affinity_sq / num_swap_candidates
    var_affinity = mean_affinity_sq - mean_affinity * mean_affinity
    var_affinity = max(var_affinity, 1e-12)
    stdev_affinity = sqrt(var_affinity)

    # Compute swap scores
    max_score = -Inf
    @inbounds for i in 1:num_swap_candidates
        (county_id, partition_id) = swaps[i]
        @inbounds affinity = scores[i]
        affinity_zscore = (affinity - mean_affinity) / stdev_affinity
        population_score = population_scores[partition_id]
        score = (affinity_zscore + population_score) / partitioner.temperature
        max_score = max(score, max_score)
        @inbounds scores[i] = score
    end

    # Convert scores to cumulative softmax sum
    scores[1] = exp(scores[1] - max_score)
    @inbounds for i in 2:length(scores)
        @inbounds scores[i] = exp(scores[i] - max_score) + scores[i-1]
    end
    prob_denom = scores[end]

    # Pick county to swap and perform swap
    while true
        rand = Random.rand(Float64)
        i = Base.Sort.searchsortedfirst(scores, rand * prob_denom)
        (county_id, partition_id) = swaps[i]
        if can_swap_county(partitioner, county_id, partition_id)
            swap_county!(partitioner, county_id, partition_id)
            break
        end
    end

end

end  # module SimulatedAnnealing
