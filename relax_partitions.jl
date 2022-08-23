#
# Relax graph partitions
#
import DataStructures
using DataStructures: DefaultDict, SortedDict
import DelimitedFiles
import Printf
import Serialization
import StatsBase
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "graph.jl"))

# Import interaction graph
println("Importing interaction graph...")
interaction_graph = Serialization.deserialize(interaction_graph_file)

# Import geography graph
println("Importing geography graph...")
geography_graph = Serialization.deserialize(geography_graph_file)

# Import county data
println("Importing county data...")
(county_data, _) = DelimitedFiles.readdlm(
    county_data_file,
    '\t',
    header=true,
)
populations = Dict{Int64, Int64}()
total_population = 0
for row in 1:size(county_data, 1)
    id = county_data[row, 1]
    population = county_data[row, 2]
    populations[id] = population
    global total_population += population
end

# Function to compute affinity of each county to each partition
function compute_partition_affinities(
    interaction_graph::Dict{Int64, Dict{Int64, Float64}},
    county_to_partition::Dict{Int64, Int64},
    )::Dict{Int64, DefaultDict{Int64, Float64}}
    partition_affinities = Dict{Int64, DefaultDict{Int64, Float64}}(
        county => DefaultDict{Int64, Float64}(0.0)
        for county in keys(interaction_graph)
    )
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
    partition_populations::DefaultDict{Int64, Int64} # part -> pop
    partition_affinities::Dict{Int64, DefaultDict{Int64, Float64}} # county -> part -> affinity
end
function PartitionData(
    interaction_graph::Dict{Int64, Dict{Int64, Float64}},
    geography_graph::Dict{Int64, Dict{Int64, Float64}},
    county_to_partition::Dict{Int64, Int64},
    )::PartitionData

    # Load initial partition
    partitions = Set{Int64}()
    partition_to_counties = DefaultDict{Int64, Set{Int64}}(Set{Int64})
    partition_populations = DefaultDict{Int64, Int64}(0)
    for (county, partition) in county_to_partition
        push!(partitions, partition)
        push!(partition_to_counties[partition], county)
        partition_populations[partition] += populations[county]
    end

    # Compute affinity between counties and partitions
    partition_affinities = compute_partition_affinities(
        interaction_graph,
        county_to_partition,
    )

    return PartitionData(
        interaction_graph,
        geography_graph,
        partitions,
        county_to_partition,
        partition_to_counties,
        partition_populations,
        partition_affinities,
    )

end

# Import partition data
println("Importing partitions...")
partition_data = DelimitedFiles.readdlm(partition_file, '\t')
county_to_partition = Dict{Int64, Int64}()
for row in 1:size(partition_data, 1)
    county = partition_data[row, 1]
    partition = partition_data[row, 2]
    county_to_partition[county] = partition
end
partition_data = PartitionData(
    interaction_graph,
    geography_graph,
    county_to_partition,
)

function transfer_county_to_partition(
    county::Int64,
    new_partition::Int64,
    partition_data::PartitionData,
    )
    interaction_graph = partition_data.interaction_graph
    partitions = partition_data.partitions
    county_to_partition = partition_data.county_to_partition
    partition_to_counties = partition_data.partition_to_counties
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities
    old_partition = county_to_partition[county]

    # Transfer county
    county_to_partition[county] = new_partition
    delete!(partition_to_counties[old_partition], county)
    push!(partition_to_counties[new_partition], county)

    # Transfer population
    partition_populations[old_partition] -= populations[county]
    partition_populations[new_partition] += populations[county]

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
    end

end

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
        sample_weights = StatsBase.Weights([
            1 / partition_affinities[neighbor][county_to_partition[neighbor]]
            for neighbor in neighbors_collect
        ])
        county = StatsBase.sample(neighbors_collect, sample_weights)

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
        sample_weights = StatsBase.Weights([
            1 / partition_affinities[county][partition]
            for county in boundary_collect
        ])
        county = StatsBase.sample(boundary_collect, sample_weights)

        # Choose partition to recieve county
        neighbor_partitions = Set{Int64}(
            county_to_partition[neighbor]
            for neighbor in keys(geography_graph[county])
        )
        delete!(neighbor_partitions, partition)
        neighbor_partitions_collect = collect(neighbor_partitions)
        sample_weights = StatsBase.Weights([
            partition_affinities[county][neighbor_partition]
            for neighbor_partition in neighbor_partitions_collect
        ])
        new_partition = StatsBase.sample(
            neighbor_partitions_collect,
            sample_weights,
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
            partition = maximum(partitions) + 1
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

# Perform rebalancing steps
println("Rebalancing partitions...")
for iter in 1:relaxation_steps

    # Aim to evenly divide population between partitions
    target_population = total_population / num_partitions
    target_population = round(Int64, target_population)

    # Randomly pick partition to adjust
    partition_populations = partition_data.partition_populations
    partitions = collect(partition_data.partitions)
    sample_weights = StatsBase.Weights([
        abs(partition_populations[partition] - target_population)
        for partition in partitions
    ])
    partition = StatsBase.sample(partitions, sample_weights)

    # Grow or shrink partition to achieve target population
    if partition_populations[partition] < target_population
        grow_partition(
            target_population,
            partition,
            partition_data,
        )
    else
        shrink_partition(
            target_population,
            partition,
            partition_data,
        )
    end

    # Split any disconnected partitions
    split_disconnected_partitions(partition_data)

end

# Print partition populations
println("Partition populations...")
for (partition, population) in partition_data.partition_populations
    @Printf.printf("Partition %d: %d\n", partition, population)
end

# Output results to file
println("Exporting partitions...")
partition_table = Array{Int64, 2}(
    undef,
    (length(partition_data.county_to_partition), 2),
)
for (row, (county, partition)) in enumerate(partition_data.county_to_partition)
    partition_table[row, 1] = county
    partition_table[row, 2] = partition
end
DelimitedFiles.writedlm(partition_file, partition_table, '\t')
