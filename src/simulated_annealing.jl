module SimulatedAnnealing

import Base.Threads
import DataStructures
import Random

import ..DataFiles
import ..DataGraphs
import ..Graph
import ..Plot

function contiguous_partition!(
    county_to_partition::Dict{UInt, UInt},
    partition_to_counties::Dict{UInt, Set{UInt}},
    partition_center_county::Dict{UInt, UInt},
    graph::Graph.WeightedGraph,
    )

    # Iterate through partitions
    for partition_id in keys(partition_to_counties)

        # BFS from partition center
        center_id = partition_center_county[partition_id]
        found_set = Set{UInt}([center_id])
        search_queue = DataStructures.Queue{UInt}()
        DataStructures.enqueue!(search_queue, center_id)
        while !isempty(search_queue)
            current = DataStructures.dequeue!(search_queue)
            for neighbor in keys(graph[current])
                if !in(neighbor, found_set)
                    push!(found_set, neighbor)
                    DataStructures.enqueue!(search_queue, neighbor)
                end
            end
        end

        # Identify disconnected counties
        partition_counties = partition_to_counties[partition_id]
        disconnected = setdiff(partition_counties, found_set)
        setdiff!(partition_counties, disconnected)

        # Reassign disconnected counties to neighboring partitions
        while !isempty(disconnected)
            for county_id in collect(disconnected)
                for neighbor in keys(graph[county_id])
                    neighbor_partition_id = county_to_partition[neighbor]
                    if neighbor == partition_id
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
    adjacency_graph::Graph.WeightedGraph
    interaction_graph::Graph.WeightedGraph
    county_to_partition::Dict{UInt, UInt}
    partition_to_counties::Dict{UInt, Set{UInt}}
    county_populations::Dict{UInt, UInt}
    partition_populations::Dict{UInt, UInt}
    swap_candidates::Dict{UInt, Dict{UInt, Float64}}
    plotter::Plot.Plotter
end

function Partitioner(
    num_partitions::UInt,
    state_ids::AbstractVector{UInt},
    interaction_personal_stdev::Float64,
    interaction_max_distance::Float64,
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
    for (partition_id, counties) in partition_to_counties
        partition_populations[partition_id] = sum(
            [county_populations[county_id] for county_id in counties]
        )
    end

    # Candidate county swaps
    swaps_candidates = Dict{UInt, Dict{UInt, Float64}}(
        id => Dict{UInt, Float64}() for id in county_ids)

    # Plotter
    plotter = Plot.Plotter(county_to_partition)

    # Construct partitioner object
    partitioner = Partitioner(
        adjacency_graph,
        interaction_graph,
        county_to_partition,
        partition_to_counties,
        county_populations,
        partition_populations,
        swaps_candidates,
        plotter,
    )

    # Find swap candidates
    @Base.Threads.threads for county_id in county_ids
        update_county_swap_candidates!(partitioner, county_id)
    end

    return partitioner

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

end  # module SimulatedAnnealing
