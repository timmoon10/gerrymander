module SimulatedAnnealing

import Base.Threads
import DataStructures
import Random

import ..DataFiles
import ..DataGraphs
import ..Graph

mutable struct Partitioner

    # Partitioner options
    temperature::Float64
    population_weight::Float64

    # County data
    adjacency_graph::Graph.WeightedGraph
    interaction_graph::Graph.WeightedGraph
    county_populations::Dict{UInt, UInt}

    # Partition data
    county_to_partition::Dict{UInt, UInt}
    partition_to_counties::Dict{UInt, Set{UInt}}
    partition_populations::Dict{UInt, UInt}
    swap_candidates::Dict{UInt, Dict{UInt, Float64}}

    # Log-interpolation of partitioner options
    interp_step::UInt
    interp_max_step::UInt
    interp_log_temperature_start::Float64
    interp_log_temperature_end::Float64
    interp_log_population_weight_start::Float64
    interp_log_population_weight_end::Float64

    # Function to parse user commands
    parse_command_func::Union{Function, Nothing}

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
        (county_to_partition, partition_to_counties) = Graph.random_partition(
            num_partitions,
            adjacency_graph,
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
        partitioner = new(
            temperature,
            population_weight,
            adjacency_graph,
            interaction_graph,
            county_populations,
            county_to_partition,
            partition_to_counties,
            partition_populations,
            swaps_candidates,
            0,
            0,
            0,
            0,
            0,
            0,
            nothing,
        )

        # Register logic for user commands
        partitioner.parse_command_func = _make_parse_command_func(partitioner)

        # Find swap candidates
        @Base.Threads.threads for county_id in county_ids
            update_county_swap_candidates!(partitioner, county_id)
        end

        return partitioner

    end

end

function _make_parse_command_func(partitioner::Partitioner)::Function

    "Logic for user commands"
    function parse_command_func(command::String)

        # Help message
        if command == "help"
            println()
            println("Commands")
            println("--------")
            println("help: help message")
            println("exit: exit")
            println("pause: pause animation")
            println("info: partitioner state")
            println("reset: reset partitioner properties")
            println("interp: start/stop property interpolation")
            println()
            println("Properties")
            println("----------")
            println("temperature")
            println("population weight")
            println("interp steps")
            println("interp temperature")
            println("interp population weight")
            println()
            return
        end

        # Partitioner info
        if command == "info"
            println()
            print_info(partitioner)
            println()
            return
        end

        # Reset partitioner
        if command == "reset"
            partitioner.temperature = 1.0
            partitioner.population_weight= 1.0
            partitioner.interp_step = 0
            partitioner.interp_max_step = 0
            partitioner.interp_log_temperature_end = 0
            partitioner.interp_log_population_weight_end = 0
            return
        end

        # Start/stop property interpolation
        if command == "interp"
            if partitioner.interp_step == 0
                start_interp!(partitioner)
            else
                partitioner.interp_step = 0
            end
            return
        end

        # Partitioner properties
        command_split = split(command, "=", limit=2)
        name = strip(command_split[1])
        value = length(command_split) > 1 ? strip(command_split[2]) : ""
        if name == "temperature"
            partitioner.temperature = parse(Float64, value)
            return
        end
        if name == "population weight"
            partitioner.population_weight = parse(Float64, value)
            return
        end
        if name == "interp steps"
            partitioner.interp_max_step = parse(UInt, value)
            return
        end
        if name == "interp temperature"
            value = parse(Float64, value)
            if value <= 0
                println("Invalid interp temperature: ", value)
            else
                partitioner.interp_log_temperature_end = log(value)
            end
            return
        end
        if name == "interp population weight"
            value = parse(Float64, value)
            if value <= 0
                println("Invalid interp population weight: ", value)
            else
                partitioner.interp_log_population_weight_end = log(value)
            end
            return
        end

        # Unrecognized command
        println("Unrecognized command: ", command)

    end

    return parse_command_func

end

function print_info(partitioner::Partitioner)

    # Print partitioner state
    println("Partitioner properties")
    println("----------------------")
    println("Temperature: ", partitioner.temperature)
    println("Population weight: ", partitioner.population_weight)
    print("\n")

    # Print interpolation state
    println("Interpolation state")
    println("-------------------")
    println("Step: ", partitioner.interp_step, " / ", partitioner.interp_max_step)
    println("Target temperature: ", exp(partitioner.interp_log_temperature_end))
    println("Target population weight: ", exp(partitioner.interp_log_population_weight_end))
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

end

function start_interp!(partitioner::Partitioner)

    # Set interpolation start points
    if partitioner.temperature <= 0
        partitioner.temperature = exp(partitioner.interp_log_temperature_end)
    end
    partitioner.interp_log_temperature_start = log(partitioner.temperature)
    if partitioner.population_weight <= 0
        partitioner.population_weight = exp(partitioner.interp_log_population_weight_end)
    end
    partitioner.interp_log_population_weight_start = log(partitioner.population_weight)

    # Interpolation step 1
    partitioner.interp_step = 1

end

function interp_step!(partitioner::Partitioner)

    # Trivial cases
    step = partitioner.interp_step
    max_step = partitioner.interp_max_step
    if step == 0
        # Interpolation is not active
        return
    elseif step >= max_step
        # Interpolation has finished
        partitioner.temperature = exp(partitioner.interp_log_temperature_end)
        partitioner.population_weight = exp(partitioner.interp_log_population_weight_end)
        partitioner.interp_step = 0
        return
    end

    "Log scale interpolation"
    function log_interp(
        start_log_val::Float64,
        end_log_val::Float64,
        frac::Float64
        )::Float64
        log_val = (end_log_val - start_log_val) * frac + start_log_val
        return exp(log_val)
    end

    # Update partitioner with interpolated properties
    frac = Float64(step - 1) / (max_step - 1)
    partitioner.temperature = log_interp(
        partitioner.interp_log_temperature_start,
        partitioner.interp_log_temperature_end,
        frac,
    )
    partitioner.population_weight = log_interp(
        partitioner.interp_log_population_weight_start,
        partitioner.interp_log_population_weight_end,
        frac,
    )
    partitioner.interp_step += 1

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
    for neighbor_id in keys(partitioner.adjacency_graph[county_id])
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
    neighbors = collect(keys(partitioner.adjacency_graph[county_id]))
    src_neighbors = Set{UInt}()
    dst_neighbors = Set{UInt}()
    sizehint!(src_neighbors, length(neighbors))
    sizehint!(dst_neighbors, length(neighbors))
    for neighbor_id in neighbors
        if neighbor_id == county_id
            continue
        end
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
    neighborhood = collect(neighborhood)
    @Base.Threads.threads for county_id in neighborhood
        update_county_swap_candidates!(partitioner, county_id)
    end

end

function step!(partitioner::Partitioner)

    # Update interpolation state
    interp_step!(partitioner)

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
    function cum_softmax_numers(
        scores::Vector{Float64},
        max_score::Float64,
        )::Vector{Float64}
        out = Vector{Float64}()
        sizehint!(out, length(scores))
        last_p::Float64 = 0
        @inbounds for (i, score) in enumerate(scores)
            p = exp(score - max_score) + last_p
            push!(out, p)
            last_p = p
        end
        return out
    end
    cum_prob_numers = cum_softmax_numers(scores, max_score)
    prob_denom = cum_prob_numers[end]

    # Try performing swap with randomly chosen county
    need_swap = true
    for _ in 1:10
        rand = Random.rand(Float64)
        i = Base.Sort.searchsortedfirst(cum_prob_numers, rand * prob_denom)
        (county_id, partition_id) = swaps[i]
        if can_swap_county(partitioner, county_id, partition_id)
            swap_county!(partitioner, county_id, partition_id)
            need_swap = false
            break
        end
    end

    # Filter invalid swaps if we pick them too many times
    if need_swap

        # Check which swaps are valid
        swap_is_valid = Vector{Bool}(undef, num_swap_candidates)
        for i in 1:num_swap_candidates
            (county_id, partition_id) = swaps[i]
            swap_is_valid[i] = can_swap_county(partitioner, county_id, partition_id)
        end
        num_valid_swaps = count(swap_is_valid)
        if num_valid_swaps == 0
            return
        end

        # Filter invalid swaps
        max_score = -Inf
        valid_swaps = Vector{Tuple{UInt, UInt}}()
        valid_scores = Vector{Float64}()
        sizehint!(valid_swaps, num_valid_swaps)
        sizehint!(valid_scores, num_valid_swaps)
        for i in 1:num_valid_swaps
            if !swap_is_valid[i]
                continue
            end
            push!(valid_swaps, swaps[i])
            push!(valid_scores, scores[i])
            max_score = max(scores[i], max_score)
        end

        # Convert scores to cumulative softmax sum
        cum_prob_numers = cum_softmax_numers(valid_scores, max_score)
        prob_denom = cum_prob_numers[end]

        # Randomly pick county and perform swap
        rand = Random.rand(Float64)
        i = Base.Sort.searchsortedfirst(cum_prob_numers, rand * prob_denom)
        (county_id, partition_id) = valid_swaps[i]
        swap_county!(partitioner, county_id, partition_id)
        need_swap = false

    end

end

end  # module SimulatedAnnealing
