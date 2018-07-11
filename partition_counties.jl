#!/usr/bin/julia
#
# Partition counties based on interaction graph and geography graph
#
using ProtoBuf
import DataStructures
using Metis
include(AbstractString(dirname(@__FILE__)) * "/common.jl")

# Import county geography graph
println("Importing county geography graph...")
geography_graph_proto = Graph()
open(geography_graph_file, "r") do f
    readproto(f, geography_graph_proto)
end
geography_graph = Set{Tuple{UInt32, UInt32}}()
for edge_proto in geography_graph_proto.edge
    geoid1 = min(edge_proto.node1, edge_proto.node2)
    geoid2 = max(edge_proto.node1, edge_proto.node2)
    push!(geography_graph, (geoid1, geoid2))
end

# Import county interaction graph
println("Importing county interaction graph...")
interaction_graph_proto = Graph()
open(interaction_graph_file, "r") do f
    readproto(f, interaction_graph_proto)
end

# Determine matrix indices corresponding to GEOIDs
num_inds = length(interaction_graph_proto.node)
geoid_to_ind = Dict{UInt32, Int32}()
ind_to_geoid = Array{UInt32, 1}(num_inds)
for ind = 1:num_inds
    geoid = interaction_graph_proto.node[ind]
    geoid_to_ind[geoid] = ind
    ind_to_geoid[ind] = geoid
end

# Import county population data
println("Importing county population data...")
(county_data, _) = readdlm(county_data_file, '\t', header=true)
pops = Array{Int32, 1}(num_inds)
for row in 1:size(county_data, 1)
    geoid = Int32(county_data[row, 1])
    if haskey(geoid_to_ind, geoid)
        ind = geoid_to_ind[geoid]
        pops[ind] = Int32(county_data[row, 2])
    end
end

# Objects for COO sparse adjacency matrix
println("Constructing COO sparse adjacency matrix...")
rows = Array{Int32, 1}()
cols = Array{Int32, 1}()
vwgts = Array{Int32, 1}()
adjwgts = Array{Float32, 1}()

# Construct edges in adjacency matrix
max_weight = 0
for edge_proto in interaction_graph_proto.edge

    # Check if counties are adjacent
    geoid1 = min(edge_proto.node1, edge_proto.node2)
    geoid2 = max(edge_proto.node1, edge_proto.node2)
    if !in((geoid1, geoid2), geography_graph)
        continue
    end
        
    # Determine edge weight
    weight = edge_proto.weight
    max_weight = max(weight, max_weight)
    
    # Add current edge to adjacency matrix
    ind1 = geoid_to_ind[geoid1]
    ind2 = geoid_to_ind[geoid2]
    push!(rows, ind1)
    push!(cols, ind2)
    push!(adjwgts, weight)
    push!(rows, ind2)
    push!(cols, ind1)
    push!(adjwgts, weight)

end
num_edges = length(adjwgts)
adjwgts_int = Array{Int32, 1}(num_edges)
for i = 1:num_edges
    adjwgts_int[i] = round(Int32, max(1e6 * adjwgts[i] / max_weight, 1))
end

# Construct CSC sparse adjacency matrix
println("Constructing CSC sparse adjacency matrix...")
adj = sparse(rows, cols, adjwgts_int)

# Partition graph
println("Partitioning graph...")
options = -ones(Cint, Metis.METIS_NOPTIONS)
options[Metis.METIS_OPTION_NCUTS] = 4
options[Metis.METIS_OPTION_CONTIG] = 1
ubvec = Array(Cfloat, 1)
ubvec[1] = balance_tolerance
objval, part = partGraphKway(adj, num_partitions,
                             adjwgt=true, vwgt=pops,
                             ubvec=ubvec, options=options)
writedlm(partition_file, [ind_to_geoid part], '\t')
