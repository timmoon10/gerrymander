module DataGraphs

import Base.Iterators
import DataStructures
import Serialization

import ..Constants
import ..DataFiles

function county_adjacency_graph(
    state_ids::AbstractVector{UInt},
    )::Dict{UInt, Dict{UInt, Float64}}

    # Return immediately if cached data exists
    cache_file = joinpath(
        DataFiles.root_dir(),
        "results",
        "county_adjacency.bin",
    )
    if isfile(cache_file)
        return Serialization.deserialize(cache_file)
    end

    # Load county data
    county_data = DataFiles.load_county_populations(state_ids)
    county_ids = Set{UInt}(county_data[:,1])

    # Load county boundaries
    county_boundaries = DataFiles.load_county_boundaries()
    filter!(p -> in(p.first, county_ids), county_boundaries)

    # Empty graph
    println("Constructing county adjacency graph...")
    graph = Dict{UInt, Dict{UInt, Float64}}(
        id => Dict{UInt, Float64}() for id in county_ids)

    # Compute contribution of border segments to geography graph
    Coord = Vector{Float64}
    Segment = Set{Coord}
    segment_owners = Dict{Segment, UInt}()
    deg_to_rad::Float64 = pi / 180
    earth_radius = Constants.earth_radius
    @inbounds for id in county_ids
        multipolygon = county_boundaries[id]
        @inbounds for polygon in multipolygon
            @inbounds for line in polygon
                @inbounds for i in 1:length(line)-1
                    segment = Segment([line[i], line[i+1]])
                    if !haskey(segment_owners, segment)
                        # Wait for another region to encounter segment
                        segment_owners[segment] = id
                    else

                        # Another region shares this border segment
                        neighbor = segment_owners[segment]
                        delete!(segment_owners, segment)
                        if neighbor == id
                            continue
                        end

                        # Compute segment length with Haversine formula
                        @inbounds l1 = deg_to_rad * line[i][1]
                        @inbounds ph1 = deg_to_rad * line[i][2]
                        @inbounds l2 = deg_to_rad * line[i+1][1]
                        @inbounds ph2 = deg_to_rad * line[i+1][2]
                        hav_angle = (
                            sin((ph2-ph1)/2)^2
                            + cos(ph1) * cos(ph2) * sin((l2-l1)/2)^2
                        )
                        d = 2 * earth_radius * asin(sqrt(hav_angle))

                        # Add segment length to graph
                        if !haskey(graph[id], neighbor)
                            graph[id][neighbor] = 0.0
                        end
                        graph[id][neighbor] += d
                        graph[neighbor][id] = graph[id][neighbor]

                    end
                end
            end
        end
    end

    # Write results to file
    Serialization.serialize(cache_file, graph)
    println("Saved county adjacency graph at $cache_file...")
    return graph

end

function county_interaction_graph(
    state_ids::AbstractVector{UInt},
    personal_stdev::Float64,
    max_distance::Float64,
    )::Dict{UInt, Dict{UInt, Float64}}

    # Return immediately if cached data exists
    cache_file = joinpath(
        DataFiles.root_dir(),
        "results",
        "county_interaction.bin",
    )
    if isfile(cache_file)
        return Serialization.deserialize(cache_file)
    end

    # Load county data
    county_data = DataFiles.load_county_populations(state_ids)
    county_ids = Set{UInt}(county_data[:,1])

    # Empty graph
    println("Constructing county interaction graph...")
    graph = Dict{UInt, Dict{UInt, Float64}}(
        id => Dict{UInt, Float64}() for id in county_ids)

    # Binned cache of county data
    CountyData = Tuple{UInt, Float64, Float64, Float64, Float64, Float64}
    binned_county_data = DataStructures.DefaultDict{
        Tuple{Int64, Int64, Int64},
        Set{CountyData},
    }(Set{CountyData})

    # Iterate through counties
    personal_var = personal_stdev ^ 2
    max_distance_sq = max_distance ^2
    @inbounds for row in 1:size(county_data, 1)

        # Current county data
        @inbounds id1::UInt = county_data[row, 1]
        @inbounds pop1::Float64 = county_data[row, 2]
        @inbounds x1::Float64 = county_data[row, 3]
        @inbounds y1::Float64 = county_data[row, 4]
        @inbounds z1::Float64 = county_data[row, 5]
        @inbounds var1::Float64 = county_data[row, 6]

        # Current bin
        bin = (
            floor(Int64, x1 / max_distance),
            floor(Int64, y1 / max_distance),
            floor(Int64, z1 / max_distance),
        )

        # Search for counties in nearby bins
        for bin in Base.Iterators.product(
            bin[1]-1:bin[1]+1,
            bin[2]-1:bin[2]+1,
            bin[3]-1:bin[3]+1,
            )
            if !haskey(binned_county_data, bin)
                continue
            end
            for (id2, pop2, x2, y2, z2, var2) in binned_county_data[bin]
                # Compute edge weight between counties
                dist2 = (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2
                if dist2 > max_distance_sq
                    continue
                end
                interaction_var = var1 + var2 + 2 * personal_var
                unit_interaction = (
                    exp(-0.5 * dist2 / interaction_var)
                    / (2 * pi * sqrt(interaction_var))
                )
                interaction_score = pop1 * pop2 * unit_interaction
                graph[id1][id2] = interaction_score
                graph[id2][id1] = interaction_score
            end
        end

        # Add current county data to binned cache
        push!(binned_county_data[bin], (id1, pop1, x1, y1, z1, var1))

    end

    # Write results to file
    Serialization.serialize(cache_file, graph)
    println("Saved county interaction graph at $cache_file...")
    return graph

end

end  # module DataGraphs
