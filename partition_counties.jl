#!/usr/bin/julia

# Parameters
project_dir = "/home/moon/Documents/gerrymander/"
output_dir  = project_dir * "/output"
props_file  = output_dir * "/props.prototxt"
graph_file  = output_dir * "/graph.prototxt"
parts_file  = output_dir * "/parts.csv"
num_partitions = 48
general_attraction = 100
partisan_attraction = 200
partisan_repulsion = 20

# Import packages
using ProtoBuf
import DataStructures
using Metis

# Read county properties and county graph from file
println("Reading county properties and graph data...")
include(output_dir * "/gerrymander_pb.jl")
props_list_proto = CountyPropertiesList()
graph_proto = Graph()
open(props_file, "r") do f
    readproto(f, props_list_proto)
end
open(graph_file, "r") do f
    readproto(f, graph_proto)
end

# Objects for COO sparse adjacency matrix
println("Constructing COO sparse adjacency matrix...")
rows = Array{Int32, 1}()
cols = Array{Int32, 1}()
vwgts = Array{Int32, 1}()
adjwgts = Array{Float32, 1}()

# Determine matrix indices corresponding to GEOIDs
num_inds = length(graph_proto.node)
geoid_to_ind = Dict{UInt32, Int32}()
ind_to_geoid = Array{UInt32, 1}(num_inds)
for ind = 1:num_inds
    geoid = graph_proto.node[ind]
    geoid_to_ind[geoid] = ind
    ind_to_geoid[ind] = geoid
end

# Read properties for each county
pops = Array{Int32, 1}(num_inds)
densities = Array{Tuple{Float32, Float32, Float32}, 1}(num_inds)
coords = Array{Float32, 2}(num_inds, 2)
for props_proto in props_list_proto.county_props
    geoid = props_proto.geoid
    if !haskey(geoid_to_ind, geoid)
        continue
    end
    ind = geoid_to_ind[geoid]
    pop_density = Float32(props_proto.population) / props_proto.area
    dem_density = Float32(props_proto.dem_votes) / props_proto.area
    gop_density = Float32(props_proto.gop_votes) / props_proto.area
    pops[ind] = props_proto.population
    densities[ind] = (pop_density, dem_density, gop_density)
    coords[ind, 1] = props_proto.interior_point.long
    coords[ind, 2] = props_proto.interior_point.lat
end

# Construct edges in adjacency matrix
max_weight = 0
for edge_proto in graph_proto.edge
    ind1 = geoid_to_ind[edge_proto.node1]
    ind2 = geoid_to_ind[edge_proto.node2]

    # Determine edge weight
    border_length = edge_proto.weight
    (pop_density1, dem_density1, gop_density1) = densities[ind1]
    (pop_density2, dem_density2, gop_density2) = densities[ind2]
    general_interaction = general_attraction * pop_density1 * pop_density2
    partisan_interaction = (partisan_attraction * (dem_density1 * dem_density2
                                                   + gop_density1 * gop_density2)
                            - partisan_repulsion * (dem_density1 * gop_density2
                                                    + gop_density1 * dem_density2))
    weight = border_length * border_length * (general_interaction + partisan_interaction)
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
ubvec[1] = 4
objval, part = partGraphKway(adj, num_partitions,
                             adjwgt=true, vwgt=pops,
                             ubvec=ubvec, options=options)
writecsv(parts_file, [ind_to_geoid part])
