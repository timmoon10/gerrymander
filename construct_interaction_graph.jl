#
# Construct county interaction matrix
#
import DataStructures
import DelimitedFiles
import Serialization
include(joinpath(dirname(@__FILE__), "common.jl"))
include(joinpath(project_dir, "graph.jl"))

# Import county population data
println("Importing county population data...")
(county_data::Array{Float64, 2}, _) = DelimitedFiles.readdlm(county_data_file, '\t', header=true)

function make_interaction_graph(
    county_data::Array{Float64, 2},
    )::Dict{Int64, Dict{Int64, Float64}}

    # Construct empty graph
    graph = Dict{Int64, Dict{Int64, Float64}}(
        Int64(id) => Dict{Int64, Float64}() for id in county_data[:,1]
    )

    # Binned cache of county data
    binned_county_data = DataStructures.DefaultDict{
        Tuple{Int64, Int64, Int64},
        Set{Array{Float64}},
    }(Set{Array{Float64}})

    personal_var::Float64 = personal_stdev ^ 2
    max_dist2::Float64 = max_dist ^ 2
    for row in 1:size(county_data, 1)

        # Current county data
        county1_data::Array{Float64} = county_data[row, :]
        id1::Int64 = county1_data[1]
        pop1 = county1_data[2]
        dem_votes1 = county1_data[3]
        gop_votes1 = county1_data[4]
        x1 = county1_data[5]
        y1 = county1_data[6]
        z1 = county1_data[7]
        var1 = county1_data[8]

        # Search for counties in nearby bins
        pos1 = [x1, y1, z1]
        bins_min = floor.(Int64, (pos1 .- max_dist) ./ bin_size)
        bins_max = floor.(Int64, (pos1 .+ max_dist) ./ bin_size)
        for bin in Iterators.product(
            bins_min[1]:bins_max[1],
            bins_min[2]:bins_max[2],
            bins_min[3]:bins_max[3],
            )
            if !haskey(binned_county_data, bin)
                continue
            end
            for county2_data in binned_county_data[bin]

                # Nearby county data
                id2::Int64 = county2_data[1]
                pop2 = county2_data[2]
                dem_votes2 = county2_data[3]
                gop_votes2 = county2_data[4]
                x2 = county2_data[5]
                y2 = county2_data[6]
                z2 = county2_data[7]
                var2 = county2_data[8]

                # Compute edge weight between counties
                dist2 = (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2
                if dist2 > max_dist2
                    continue
                end
                interaction_var = var1 + var2 + 2 * personal_var
                unit_interaction = (
                    exp(-0.5 * dist2 / interaction_var)
                    / (2 * pi * sqrt(interaction_var))
                )
                eff_pop_product = (
                    pop1 * pop2
                    + partisan_attraction
                    * (dem_votes1 * dem_votes2 + gop_votes1 * gop_votes2)
                    - partisan_repulsion
                    * (dem_votes1 * gop_votes2 + gop_votes1 * dem_votes2)
                )
                interaction_score = eff_pop_product * unit_interaction
                graph[id1][id2] = interaction_score
                graph[id2][id1] = interaction_score

            end
        end

        # Add current county data to binned cache
        bin = (
            floor(Int64, x1 / bin_size),
            floor(Int64, y1 / bin_size),
            floor(Int64, z1 / bin_size),
        )
        push!(binned_county_data[bin], county1_data)

    end

    return graph
end

# Construct interaction graph
println("Constructing interaction graph...")
graph = make_interaction_graph(county_data)

# Remove unconnected components of graph
println("Removing unconnected components of interaction graph...")
search_start = minimum(keys(graph))
connected_ids = find_connected_vertices(graph, search_start)
if length(connected_ids) != length(graph)
    println("Interaction graph is unconnected!")
    graph = construct_subgraph(graph, connected_ids)
end

# Write results to file
println("Exporting interaction graph...")
Serialization.serialize(interaction_graph_file, graph)
