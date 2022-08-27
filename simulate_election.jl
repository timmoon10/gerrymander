#
# Simulate election results
#
import DataStructures
import DelimitedFiles
import Printf
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "libgeos_utils.jl"))

# Import county data
println("Importing county data...")
(county_data, _) = DelimitedFiles.readdlm(
    county_data_file,
    '\t',
    header=true,
)
county_populations = Dict{Int64, Int64}()
county_dem_votes = Dict{Int64, Int64}()
county_gop_votes = Dict{Int64, Int64}()
for row in 1:size(county_data, 1)
    county = county_data[row, 1]
    county_populations[county] = county_data[row, 2]
    county_dem_votes[county] = county_data[row, 3]
    county_gop_votes[county] = county_data[row, 4]
end

# Import partition data
println("Importing partition data...")
partition_data = DelimitedFiles.readdlm(partition_file, '\t')
partitions = Dict{Int64, Int64}()
partition_populations = DataStructures.DefaultDict{Int64, Int64}(0)
for row in 1:size(partition_data, 1)
    county = partition_data[row, 1]
    partition = partition_data[row, 2]
    partitions[county] = partition
    partition_populations[partition] += county_populations[county]
end

# Huntington-Hill apportionment method for US House
function house_apportionment(state_populations)
    # House seats, excluding Hawaii and Alaska
    seats_remaining = 435 - 3
    state_seats = Dict{Int64, Int64}()
    priorities = DataStructures.BinaryMaxHeap{Tuple{Float64,Int64}}()
    for (state, population) in state_populations
        state_seats[state] = 1
        push!(priorities, (population / sqrt(2), state))
    end

    # House seats remaining, excluding Hawaii and Alaska
    seats_to_allocate = 435 - 3 - length(state_populations)
    for _ in 1:seats_to_allocate
        (_, state) = pop!(priorities)
        seats = state_seats[state] + 1
        state_seats[state] = seats
        priority = state_populations[state] / sqrt(seats * (seats+1))
        push!(priorities, (priority, state))
    end

    return state_seats
end

# Compute electoral college votes
partition_electoral_votes = house_apportionment(partition_populations)
for partition in keys(partition_electoral_votes)
    partition_electoral_votes[partition] += 2
end

# Compute per-partition data
partition_populations = DataStructures.DefaultDict{Int64, Int64}(0)
partition_dem_votes = DataStructures.DefaultDict{Int64, Int64}(0)
partition_gop_votes = DataStructures.DefaultDict{Int64, Int64}(0)
for (county, partition) in partitions
    partition_populations[partition] += county_populations[county]
    partition_dem_votes[partition] += county_dem_votes[county]
    partition_gop_votes[partition] += county_gop_votes[county]
end

# Print results
electoral_votes = zeros(Int64, 2)
partition_leans = zeros(Int64, 4)
for partition in keys(partition_populations)
    dem_vote = partition_dem_votes[partition]
    gop_vote = partition_gop_votes[partition]
    total_vote = dem_vote + gop_vote
    lean = ""
    if dem_vote < 0.45 * total_vote
        lean = "R"
        partition_leans[1] += 1
    elseif dem_vote < 0.5 * total_vote
        lean = "Lean R"
        partition_leans[2] += 1
    elseif dem_vote <= 0.55 * total_vote
        lean = "Lean D"
        partition_leans[3] += 1
    else
        lean = "D"
        partition_leans[4] += 1
    end
    if dem_vote < gop_vote
        electoral_votes[1] += partition_electoral_votes[partition]
    else
        electoral_votes[2] += partition_electoral_votes[partition]
    end
    Printf.@printf(
        "Partition %d (%s): population = %d, electoral votes = %d, DEM = %d, GOP = %d\n",
        partition,
        lean,
        partition_populations[partition],
        partition_electoral_votes[partition],
        dem_vote,
        gop_vote,
    )
end
println("Solid R: ",partition_leans[1])
println("Lean R: ",partition_leans[2])
println("Lean D: ",partition_leans[3])
println("Solid D: ",partition_leans[4])
println("Electoral R: ",electoral_votes[1])
println("Electoral D: ",electoral_votes[2])
