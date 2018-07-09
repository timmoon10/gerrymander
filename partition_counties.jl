#!/usr/bin/julia
using ProtoBuf
import DataStructures
using Metis

# Parameters
project_dir      = chomp(readall(`git rev-parse --show-toplevel`)) * "/"
protoc           = "/home/moon/src/protobuf/src/protoc"
results_dir      = project_dir * "/results/"
county_data_file = results_dir * "/county_data.tsv"
graph_file       = results_dir * "/graph.prototxt"
parts_file       = results_dir * "/parts.tsv"
num_partitions   = 6

# Initialize protobuf
println("Initializing protobuf...")
isdir(results_dir) || mkdir(results_dir)
if !isfile(results_dir * "/gerrymander_pb.jl")
    run(`$protoc --plugin=$julia_protobuf_dir/plugin/protoc-gen-julia 
         -I=$project_dir --julia_out=$results_dir
         gerrymander.proto`)
end
include(results_dir * "/gerrymander_pb.jl")

# Import county graph
println("Importing county interaction graph...")
graph_proto = Graph()
open(graph_file, "r") do f
    readproto(f, graph_proto)
end

# Determine matrix indices corresponding to GEOIDs
num_inds = length(graph_proto.node)
geoid_to_ind = Dict{UInt32, Int32}()
ind_to_geoid = Array{UInt32, 1}(num_inds)
for ind = 1:num_inds
    geoid = graph_proto.node[ind]
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
for edge_proto in graph_proto.edge
    ind1 = geoid_to_ind[edge_proto.node1]
    ind2 = geoid_to_ind[edge_proto.node2]

    # Determine edge weight
    weight = edge_proto.weight
    max_weight = max(weight, max_weight)
    
    # Add current edge to adjacency matrix
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
    adjwgts_int[i] = max(1, round(Int32, 1e6 * adjwgts[i] / max_weight))
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
ubvec[1] = 1.25
objval, part = partGraphKway(adj, num_partitions,
                             adjwgt=true, vwgt=pops,
                             ubvec=ubvec, options=options)
writedlm(parts_file, [ind_to_geoid part], '\t')
