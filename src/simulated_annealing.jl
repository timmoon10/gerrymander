module SimulatedAnnealing

import Base.Threads
import DataStructures
import DelimitedFiles
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
        county_ids = Vector{UInt}(county_population_data[:, 1])
        county_populations = Dict{UInt, UInt}(
            county_ids[i] => county_population_data[i, 2]
            for i in 1:size(county_population_data, 1))

        # Initial partition
        (county_to_partition, partition_to_counties) = Graph.random_partition(
            num_partitions,
            adjacency_graph,
        )

        # Construct partitioner object
        partitioner = new(
            1.0,
            1.0,
            adjacency_graph,
            interaction_graph,
            county_populations,
            county_to_partition,
            partition_to_counties,
            Dict{UInt, UInt}(),
            Dict{UInt, Dict{UInt, Float64}}(),
            0,
            0,
            0,
            0,
            0,
            0,
            nothing,
        )

        # Reset partitioner data
        reset_partitions!(partitioner)

        # Register logic for user commands
        partitioner.parse_command_func = _make_parse_command_func(partitioner)

        return partitioner

    end

end

function reset_partitions!(partitioner::Partitioner)

    # Get counties and partitions
    county_to_partition = partitioner.county_to_partition
    county_ids = Vector{UInt}()
    sizehint!(county_ids, length(county_to_partition))
    partition_ids = Set{UInt}()
    for (county_id, partition_id) in county_to_partition
        push!(county_ids, county_id)
        push!(partition_ids, partition_id)
    end

    # Make sure partition data is consistent
    partition_to_counties = partitioner.partition_to_counties
    empty!(partition_to_counties)
    sizehint!(partition_to_counties, length(partition_ids))
    for partition_id in partition_ids
        partition_to_counties[partition_id] = Set{UInt}()
    end
    for (county_id, partition_id) in county_to_partition
        push!(partition_to_counties[partition_id], county_id)
    end

    # Make sure partitions are contiguous
    partition_centers = Dict{UInt, UInt}(
        id => first(partition_to_counties[id]) for id in partition_ids)
    Graph.make_partition_contiguous!(
        county_to_partition,
        partition_to_counties,
        partition_centers,
        partitioner.adjacency_graph,
    )

    # Compute partition populations
    partition_populations = partitioner.partition_populations
    county_populations = partitioner.county_populations
    empty!(partition_populations)
    sizehint!(partition_populations, length(partition_ids))
    for partition_id in partition_ids
        partition_populations[partition_id] = 0
    end
    for (county_id, partition_id) in county_to_partition
        partition_populations[partition_id] += county_populations[county_id]
    end

    # Compute swap candidates
    swap_candidates = partitioner.swap_candidates
    empty!(swap_candidates)
    sizehint!(swap_candidates, length(county_ids))
    for county_id in county_ids
        swap_candidates[county_id] = Dict{UInt, Float64}()
    end
    @Base.Threads.threads for county_id in county_ids
        update_county_swap_candidates!(partitioner, county_id)
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
            println("save: save partitions to file")
            println("load: load partitions from file")
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
            reset_partitions!(partitioner)
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

        # Parametrized commands
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
        if name == "save"
            save_partition(partitioner, String(value))
            return
        end
        if name == "load"
            load_partition!(partitioner, String(value))
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
                if in(neighbor_id, partition_counties)
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

    # Number of candidate swaps
    num_swap_candidates::Int = 0
    for county_swaps in values(partitioner.swap_candidates)
        num_swap_candidates += length(county_swaps)
    end

    # Flatten candidate swaps and compute sums
    swaps = Vector{Tuple{UInt, UInt}}()
    affinities = Vector{Float64}()
    pop_diffs = Vector{Float64}()
    sizehint!(swaps, num_swap_candidates)
    sizehint!(affinities, num_swap_candidates)
    sizehint!(pop_diffs, num_swap_candidates)
    sum_affinity::Float64 = 0
    sum_affinity_sq::Float64 = 0
    sum_pop_diff::Float64 = 0
    sum_pop_diff_sq::Float64 = 0
    for (county_id, county_swaps) in partitioner.swap_candidates
        for (dst_partition, affinity) in county_swaps
            src_partition = partitioner.county_to_partition[county_id]
            pop_diff = (
                partitioner.partition_populations[src_partition]
                - partitioner.partition_populations[dst_partition]
            )
            push!(swaps, (county_id, dst_partition))
            push!(affinities, affinity)
            push!(pop_diffs, pop_diff)
            sum_affinity += affinity
            sum_affinity_sq += affinity * affinity
            sum_pop_diff += pop_diff
            sum_pop_diff_sq += pop_diff * pop_diff
        end
    end
    num_swap_candidates = length(swaps)
    if num_swap_candidates == 0
        return
    end

    # Compute statistics
    function compute_statistics(
        sum::Float64,
        sum_sq::Float64,
        count::Int,
        )::Tuple{Float64, Float64}
        mean = sum / count
        mean_sq = sum_sq / count
        var = mean_sq - mean * mean
        var = max(var, 1e-12)
        stdev = sqrt(var)
        return (mean, stdev)
    end
    function zscore(val::Float64, stats::Tuple{Float64, Float64})::Float64
        (mean, stdev) = stats
        return (val - mean) / stdev
    end
    affinity_stats = compute_statistics(
        sum_affinity,
        sum_affinity_sq,
        num_swap_candidates,
    )
    pop_diff_stats = compute_statistics(
        sum_pop_diff,
        sum_pop_diff_sq,
        num_swap_candidates,
    )

    # Compute scores
    scores = Vector{Float64}()
    sizehint!(scores, num_swap_candidates)
    max_score = -Inf
    @inbounds for i in 1:num_swap_candidates
        @inbounds affinity = affinities[i]
        @inbounds pop_diff = pop_diffs[i]
        score = (
            zscore(affinity, affinity_stats)
            + partitioner.population_weight * zscore(pop_diff, pop_diff_stats)
        )
        score /= partitioner.temperature
        max_score = max(score, max_score)
        push!(scores, score)
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
        swap_is_valid = zeros(Bool, num_swap_candidates)
        for i in 1:num_swap_candidates
            (county_id, partition_id) = swaps[i]
            swap_is_valid[i] = can_swap_county(partitioner, county_id, partition_id)
        end
        num_valid_swaps = count(swap_is_valid)
        if num_valid_swaps == 0
            return
        end

        # Filter invalid swaps
        valid_swaps = Vector{Tuple{UInt, UInt}}()
        valid_scores = Vector{Float64}()
        sizehint!(valid_swaps, num_valid_swaps)
        sizehint!(valid_scores, num_valid_swaps)
        max_score = -Inf
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

function save_partition(
    partitioner::Partitioner,
    file::String,
    )
    county_to_partition = partitioner.county_to_partition
    county_ids = collect(keys(county_to_partition))
    sort!(county_ids)
    data = Array{UInt, 2}(undef, length(county_ids), 2)
    for (i, county_id) in enumerate(county_ids)
        data[i, 1] = county_id
        data[i, 2] = county_to_partition[county_id]
    end
    DelimitedFiles.writedlm(file, data, ',')
end

function load_partition!(
    partitioner::Partitioner,
    file::String,
    )
    data = DelimitedFiles.readdlm(file, ',', UInt)
    for i in 1:size(data, 1)
        county_id = data[i, 1]
        partition_id = data[i, 2]
        partitioner.county_to_partition[county_id] = partition_id
    end
    reset_partitions!(partitioner)
end

end  # module SimulatedAnnealing
