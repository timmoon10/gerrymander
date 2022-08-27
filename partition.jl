#
# Utility functions for manipulating partitions
#
using DataStructures: DefaultDict, SortedDict
import LinearAlgebra
import Statistics
import StatsBase
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "graph.jl"))
include(joinpath(project_dir, "math.jl"))

# Function to compute affinity of each county to each partition
function update_partition_affinities(
    partition_affinities::Dict{Int64, DefaultDict{Int64, Float64}},
    interaction_graph::Dict{Int64, Dict{Int64, Float64}},
    county_to_partition::Dict{Int64, Int64},
    )::Dict{Int64, DefaultDict{Int64, Float64}}
    empty!(partition_affinities)
    for county in keys(interaction_graph)
        partition_affinities[county] = DefaultDict{Int64, Float64}(0.0)
    end
    for (county, neighborhood) in interaction_graph
        for (neighbor, affinity) in neighborhood
            partition = county_to_partition[neighbor]
            partition_affinities[county][partition] += affinity
        end
    end
    return partition_affinities
end

# Helper struct to store partition assignments
struct PartitionData
    interaction_graph::Dict{Int64, Dict{Int64, Float64}}
    geography_graph::Dict{Int64, Dict{Int64, Float64}}
    partitions::Set{Int64}
    county_to_partition::Dict{Int64, Int64}
    partition_to_counties::DefaultDict{Int64, Set{Int64}}
    county_populations::Dict{Int64, Int64}
    partition_populations::DefaultDict{Int64, Int64} # part -> pop
    partition_affinities::Dict{Int64, DefaultDict{Int64, Float64}} # county -> part -> affinity
end
function PartitionData(
    interaction_graph::Dict{Int64, Dict{Int64, Float64}},
    geography_graph::Dict{Int64, Dict{Int64, Float64}},
    county_to_partition::Dict{Int64, Int64},
    county_populations::Dict{Int64, Int64},
    )::PartitionData
    partition_data = PartitionData(
        interaction_graph,
        geography_graph,
        Set{Int64}(), # partitions
        Dict{Int64, Int64}(), # county_to_partition
        DefaultDict{Int64, Set{Int64}}(Set{Int64}), # partition_to_counties
        county_populations,
        DefaultDict{Int64, Int64}(0), # partition_populations
        Dict{Int64, DefaultDict{Int64, Float64}}(), # partition_affinities
    )
    update_partition_data(county_to_partition, partition_data)
    return partition_data
end

function update_partition_data(
    county_to_partition::Dict{Int64, Int64},
    partition_data::PartitionData,
    )
    partitions = partition_data.partitions
    partition_to_counties = partition_data.partition_to_counties
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities

    # Clear old data
    empty!(partitions)
    empty!(partition_data.county_to_partition)
    empty!(partition_to_counties)
    empty!(partition_populations)
    empty!(partition_affinities)

    # Initialize partitions
    for (county, partition) in county_to_partition
        partition_data.county_to_partition[county] = partition
        push!(partitions, partition)
        push!(partition_to_counties[partition], county)
        partition_populations[partition] += county_populations[county]
    end

    # Compute affinity between counties and partitions
    update_partition_affinities(
        partition_affinities,
        interaction_graph,
        county_to_partition,
    )

end

function transfer_county_to_partition(
    county::Int64,
    new_partition::Int64,
    partition_data::PartitionData,
    )
    interaction_graph = partition_data.interaction_graph
    partitions = partition_data.partitions
    county_to_partition = partition_data.county_to_partition
    partition_to_counties = partition_data.partition_to_counties
    county_populations = partition_data.county_populations
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities
    old_partition = county_to_partition[county]

    # Transfer county
    push!(partitions, new_partition)
    county_to_partition[county] = new_partition
    delete!(partition_to_counties[old_partition], county)
    push!(partition_to_counties[new_partition], county)

    # Transfer population
    partition_populations[old_partition] -= county_populations[county]
    partition_populations[new_partition] += county_populations[county]

    # Transfer affinity
    for (neighbor, affinity) in interaction_graph[county]
        partition_affinities[neighbor][old_partition] -= affinity
        partition_affinities[neighbor][new_partition] += affinity
    end

    # Remove old partition if empty
    if isempty(partition_to_counties[old_partition])
        delete!(partitions, old_partition)
        delete!(partition_to_counties, old_partition)
        delete!(partition_populations, old_partition)
    elseif partition_populations[old_partition] == 0
        error("Found a partition with zero population")
    end

end

# Take weakly-held counties from neighbors until target population is
# reached
function grow_partition(
    target_population::Int64,
    partition::Int64,
    partition_data::PartitionData,
    )
    geography_graph = partition_data.geography_graph
    county_to_partition = partition_data.county_to_partition
    partition_to_counties = partition_data.partition_to_counties
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities

    # Find counties adjacent to partition
    neighbors = Set{Int64}()
    for county in partition_to_counties[partition]
        for neighbor in keys(geography_graph[county])
            if !in(neighbor, partition_to_counties[partition])
                push!(neighbors, neighbor)
            end
        end
    end

    # Add counties until partition size reaches target
    while partition_populations[partition] < target_population

        # Choose county to take
        neighbors_collect = collect(neighbors)
        weights = [
            (
                partition_affinities[county][partition]
                - partition_affinities[county][county_to_partition[county]]
            )
            for county in neighbors_collect
        ]
        zscore!(weights)
        weights = StatsBase.Weights(exp.(weights))
        county = StatsBase.sample(neighbors_collect, weights)

        # Transfer county to partition
        transfer_county_to_partition(
            county,
            partition,
            partition_data,
        )

        # Update list of adjacent counties
        delete!(neighbors, county)
        for neighbor in keys(geography_graph[county])
            if !in(neighbor, partition_to_counties[partition])
                push!(neighbors, neighbor)
            end
        end

    end

end

# Give weakly-held counties to neighbors until target population is
# reached
function shrink_partition(
    target_population::Int64,
    partition::Int64,
    partition_data::PartitionData,
    )
    geography_graph = partition_data.geography_graph
    partitions = partition_data.partitions
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities
    county_to_partition = partition_data.county_to_partition
    partition_counties = partition_data.partition_to_counties[partition]

    # Find non-interior counties
    boundary = Set{Int64}()
    for county in partition_counties
        if any(!in(neighbor, partition_counties)
               for neighbor in keys(geography_graph[county]))
            push!(boundary, county)
        end
    end

    # Remove counties until partition size reaches target
    while partition_populations[partition] > target_population

        # Choose county to eject
        boundary_collect = collect(boundary)
        weights = [
            (
                sum(aff for (_, aff) in partition_affinities[county])
                - 2*partition_affinities[county][partition]
            )
            for county in boundary_collect
        ]
        zscore!(weights)
        weights = StatsBase.Weights(exp.(weights))
        county = StatsBase.sample(boundary_collect, weights)

        # Choose partition to recieve county
        neighbor_partitions = Set{Int64}(
            county_to_partition[neighbor]
            for neighbor in keys(geography_graph[county])
        )
        delete!(neighbor_partitions, partition)
        neighbor_partitions_collect = collect(neighbor_partitions)
        weights = [
            partition_affinities[county][neighbor_partition]
            for neighbor_partition in neighbor_partitions_collect
        ]
        zscore!(weights)
        weights = StatsBase.Weights(exp.(weights))
        new_partition = StatsBase.sample(
            neighbor_partitions_collect,
            weights,
        )

        # Transfer county to other partition
        transfer_county_to_partition(
            county,
            new_partition,
            partition_data,
        )

        # Exit immediately if partition is empty
        if !in(partition, partitions)
            return
        end

        # Update list of non-interior counties
        delete!(boundary, county)
        for neighbor in keys(geography_graph[county])
            if in(neighbor, partition_counties)
                push!(boundary, neighbor)
            end
        end

    end

end

# Divide disconnected partitions into separate partitions
function split_disconnected_partitions(
    partition_data::PartitionData,
    )
    partitions = partition_data.partitions
    partition_to_counties = partition_data.partition_to_counties
    geography_graph = partition_data.geography_graph

    for partition in collect(partitions)
        partition_counties = copy(partition_to_counties[partition])
        partition_graph = construct_subgraph(
            geography_graph,
            partition_counties,
        )
        connected_counties = find_connected_vertices(
            partition_graph,
            first(partition_counties),
        )
        while length(connected_counties) != length(partition_counties)
            partition = 1
            while in(partition, partitions)
                partition += 1
            end
            push!(partitions, partition)
            setdiff!(partition_counties, connected_counties)
            for county in partition_counties
                transfer_county_to_partition(
                    county,
                    partition,
                    partition_data,
                )
            end
            partition_graph = construct_subgraph(
                partition_graph,
                partition_counties,
            )
            connected_counties = find_connected_vertices(
                partition_graph,
                first(partition_counties),
            )
        end
    end

end

# Create a new partition from within and grow
function schism_partition(
    partition::Int64,
    partition_data::PartitionData,
    )
    partitions = partition_data.partitions
    partition_counties = partition_data.partition_to_counties[partition]
    population = partition_data.partition_populations[partition]
    partition_affinities = partition_data.partition_affinities
    new_partition = maximum(partitions) + 1
    target_population = round(
        Int64,
        (rand()/2+1/4) * population,
    )

    # Choose county to eject
    counties_collect = collect(partition_counties)
    weights = [
        partition_affinities[county][partition]
        for county in counties_collect
    ]
    zscore!(weights)
    weights = StatsBase.Weights(exp.(weights))
    county = StatsBase.sample(counties_collect, weights)

    # Grow partition from ejected county
    transfer_county_to_partition(
        county,
        new_partition,
        partition_data,
    )
    grow_partition(
        target_population,
        new_partition,
        partition_data,
    )

end

function softmax_relaxation(
    target_population::Float64,
    steps::Int64,
    partition_data::PartitionData,
    )
    interaction_graph = partition_data.interaction_graph
    num_counties = length(partition_data.county_to_partition)
    local num_partitions = length(partition_data.partitions)

    # Convert county and partition IDs to local indices
    col_to_county = collect(keys(partition_data.county_to_partition))
    row_to_partition = collect(partition_data.partitions)
    county_to_col = Dict{Int64, Int64}(
        county => col
        for (col, county) in enumerate(col_to_county)
    )
    partition_to_row = Dict{Int64, Int64}(
        partition => row
        for (row, partition) in enumerate(row_to_partition)
    )

    # Construct vector with county populations
    total_population = 0
    county_populations = Array{Float64}(undef, num_counties)
    for (col, county) in enumerate(col_to_county)
        pop = partition_data.county_populations[county]
        total_population += pop
        county_populations[col] = pop
    end

    # Normalize interaction graph
    normalized_graph = Array{Dict{Int64,Float64}}(undef, num_counties)
    for col in 1:num_counties
        county = col_to_county[col]
        neighbors = collect(keys(interaction_graph[county]))
        affinities = [
            interaction_graph[county][neighbor]
            for neighbor in neighbors
        ]
        zscore!(affinities)
        normalized_graph[col] = Dict{Int64,Float64}(
            county_to_col[neighbor] => affinity
            for (neighbor, affinity) in zip(neighbors, affinities)
        )
    end

    # Initialize partition membership as one-hot
    loyalties = zeros(num_partitions, num_counties)
    for (county, partition) in partition_data.county_to_partition
        row = partition_to_row[partition]
        col = county_to_col[county]
        loyalties[row, col] = 1.0
    end
    new_loyalties = Array{Float64, 2}(undef, num_partitions, num_counties)

    for step in 1:steps

        # Compute factors for rebalancing partition populations
        partition_populations = loyalties * county_populations
        scales::Array{Float64,1} = [
            min(1.0, target_population / pop)
            for pop in partition_populations
        ]
        one_minus_scales::Array{Float64,1} = 1.0 .- scales
        grow_factors::Array{Float64,1} = [
            max(0.0, target_population - pop)
            for pop in partition_populations
        ]
        grow_factors ./= max(sum(grow_factors), 1e-8)

        # Rebalance partition populations
        for col in 1:num_counties
            loyalties_col = view(loyalties, :, col)
            transfer_population = LinearAlgebra.dot(
                one_minus_scales,
                loyalties_col)
            loyalties_col .*= scales
            loyalties_col .+= transfer_population .* grow_factors
            for row in 1:num_partitions
                loyalties_col[row] += 1e-2*randn()
            end
        end

        # Recompute partition membership
        new_loyalties = zeros(size(loyalties))
        for col in 1:num_counties
            new_loyalties_col = view(new_loyalties, :, col)
            for (neighbor_col, affinity) in normalized_graph[col]
                new_loyalties_col .+= affinity .* @view loyalties[:,neighbor_col]
            end
            softmax!(new_loyalties_col)
        end
        loyalties, new_loyalties = new_loyalties, loyalties

    end

    # Choose partition membership based on greatest loyalty
    county_to_partition = argmax(loyalties, dims=1)
    county_to_partition = Dict{Int64, Int64}(
        col_to_county[inds[2]] => row_to_partition[inds[1]]
        for inds in county_to_partition
    )

    # Update partition data
    update_partition_data(county_to_partition, partition_data)

end
