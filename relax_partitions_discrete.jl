#
# Relax graph partitions
#
import DelimitedFiles
import Printf
import Serialization
import StatsBase
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "partition.jl"))

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
county_populations = Dict{Int64, Int64}()
total_population = 0
for row in 1:size(county_data, 1)
    id = county_data[row, 1]
    population = county_data[row, 2]
    county_populations[id] = population
    global total_population += population
end

# Initialize partitions
county_to_partition = Dict{Int64, Int64}()
if isfile(partition_file)
    println("Initializing with imported partitions...")
    partition_data = DelimitedFiles.readdlm(partition_file, '\t')
    for row in 1:size(partition_data, 1)
        county = partition_data[row, 1]
        partition = partition_data[row, 2]
        county_to_partition[county] = partition
    end
else
    println("Initializing with trivial partition...")
    for county in keys(geography_graph)
        county_to_partition[county] = 1
    end
end
partition_data = PartitionData(
    interaction_graph,
    geography_graph,
    county_to_partition,
    county_populations,
)

# Print parameters
@Printf.printf("Counties: %d\n", length(interaction_graph))
@Printf.printf("Total population: %d\n", total_population)
@Printf.printf("Target partitions: %d\n", num_partitions)
@Printf.printf(
    "Target population: %d\n",
    round(Int64, total_population / num_partitions),
)
@Printf.printf("Relaxation steps: %d\n", relaxation_steps)

# Perform rebalancing steps
println("Rebalancing partitions...")
for iter in 1:relaxation_steps
    partitions = partition_data.partitions
    partition_to_counties = partition_data.partition_to_counties
    partition_populations = partition_data.partition_populations
    partition_affinities = partition_data.partition_affinities

    # Update affinities
    update_partition_affinities(
        partition_affinities,
        partition_data.interaction_graph,
        partition_data.county_to_partition,
    )

    # Aim to evenly divide population between partitions
    target_population = total_population / num_partitions
    target_population = round(Int64, target_population)

    # Randomly pick partition to adjust
    partitions_collect = collect(partitions)
    sample_weights = StatsBase.Weights([
        abs(partition_populations[partition] - target_population)
        for partition in partitions_collect
    ])
    partition = StatsBase.sample(partitions_collect, sample_weights)
    population = partition_populations[partition]

    # Manipulate partitions to achieve target population
    if population <= target_population
        # Grow partition if somewhat small
        ratio = population / target_population
        target_population = round(
            Int64,
            population + rand() * (target_population - population),
        )
        grow_partition(target_population, partition, partition_data)
    else
        ratio = target_population / population
        if length(partition_to_counties[partition]) > 1
            if (length(partitions) == 1
                || (ratio < 0.5 && rand() > 2*ratio))
                # Schism partition if very large
                schism_partition(partition, partition_data)
            else
                # Shrink partition if somewhat large
                target_population = round(
                    Int64,
                    target_population + rand() * (population - target_population),
                )
                shrink_partition(
                    target_population,
                    partition,
                    partition_data,
                )
            end
        end
    end

    # Split any disconnected partitions
    split_disconnected_partitions(partition_data)

    # Randomly destroy small partitions
    target_population = total_population / num_partitions
    target_population = round(Int64, target_population)
    for partition in collect(partitions)
        population = partition_populations[partition]
        ratio = population / target_population
        if rand() > (2*ratio)^2
            shrink_partition(0, partition, partition_data)
        end
    end

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
