module Election

using DataStructures

using ..DataFiles

"Huntington-Hill apportionment method for US House"
function house_apportionment(state_populations::Dict{UInt, UInt})

    # Allocate one seat to each state
    state_seats = Dict{UInt, UInt}(
        state => 1 for state in keys(state_populations))

    # Seat priority for each state
    priorities = DataStructures.BinaryMaxHeap{Tuple{Float64,UInt}}()
    sizehint!(priorities, length(state_populations))
    for (state, population) in state_populations
        push!(priorities, (Float64(population) / sqrt(2), state))
    end

    # Allocate remaining seats, excluding Hawaii and Alaska
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

function simulate_election(
    partition_ids::Vector{UInt},
    county_to_partition::Dict{UInt, UInt},
    partition_populations::Dict{UInt, UInt},
    )

    # Load election data
    election_data = DataFiles.load_election_data()

    # Compute votes in each partition
    dem_votes = Dict{UInt, UInt}(id => 0 for id in partition_ids)
    gop_votes = Dict{UInt, UInt}(id => 0 for id in partition_ids)
    for (county, partition) in county_to_partition
        dem_votes[partition] += election_data[county][1]
        gop_votes[partition] += election_data[county][2]
    end

    # Compute GOP fraction
    gop_fraction = Dict{UInt, Float64}()
    sizehint!(gop_fraction, length(partition_ids))
    for partition in partition_ids
        gop_fraction[partition] = (
            Float64(gop_votes[partition])
            / (dem_votes[partition] + gop_votes[partition])
        )
    end

    # Compute Electoral College votes
    dem_electoral_votes = 4  # Hawaii
    gop_electoral_votes = 3  # Alaska
    house_seats = house_apportionment(partition_populations)
    for partition in partition_ids
        electoral_votes = house_seats[partition] + 2
        if dem_votes[partition] > gop_votes[partition]
            dem_electoral_votes += electoral_votes
        else
            gop_electoral_votes += electoral_votes
        end
    end

    return gop_fraction, (dem_electoral_votes, gop_electoral_votes)

end

end  # module Election
