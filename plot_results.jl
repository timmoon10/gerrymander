#!/usr/bin/julia
#
# Plot partitions
#
import ProtoBuf, LibGEOS, GeoInterface, PyPlot
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Generate color scheme
println("Generating plot colors...")
color_list = [(86,180,233), (213,94,0), (0,158,115), (240,228,66),
                 (0,114,178), (204,121,167), (230,159,0)]
color_list = [(r/256, g/256, b/256) for (r,g,b) in color_list]
if length(color_list) < num_partitions
    color_list = [color_list;
                  rand(color_list, num_partitions - length(color_list))]
end

# Import partition data
println("Importing partition data...")
partition_data = readdlm(partition_file, '\t')
num_regions = size(partition_data, 1)
partition_list = Dict{Int64, Int64}()
for row in 1:num_regions
    geoid = Int64(partition_data[row, 1])
    partition = Int64(partition_data[row, 2])
    partition_list[geoid] = partition
end

# Read boundary data from protobuf
println("Importing geography data...")
region_list_proto = MultiPolygonList()
open(geography_data_file, "r") do f
    ProtoBuf.readproto(f, region_list_proto)
end

# Initialize Lambert conformal conic projection
println("Initializing map projection...")
min_long = 0.0
max_long = -180.0
min_lat = 90.0
max_lat = -90.0
for region_proto in region_list_proto.multi_polygon
    id = region_proto.id
    if haskey(partition_list, id)
        for polygon_proto in region_proto.polygon
            for long in polygon_proto.exterior_border.x
                min_long = min(long, min_long)
                max_long = max(long, max_long)
            end
            for lat in polygon_proto.exterior_border.y
                min_lat = min(lat, min_lat)
                max_lat = max(lat, max_lat)
            end
        end
    end
end
min_long *= pi / 180
max_long *= pi / 180
min_lat  *= pi / 180
max_lat  *= pi / 180
const ref_long = (min_long + max_long) / 2
const ref_lat = (min_lat + max_lat) / 2
const lambert_n = (log(cos(max_lat) * sec(min_lat))
                   / log(tan(pi/4 + min_lat/2) * cot(pi/4 + max_lat/2)))
const lambert_F = (cos(max_lat) * tan(pi/4+max_lat/2)^lambert_n)/ lambert_n
const lambert_rho_ref = lambert_F * cot(pi/4 + ref_lat/2)^lambert_n
function lambert_projection(long_deg, lat_deg)
    long = pi / 180 * long_deg
    lat  = pi / 180 * lat_deg
    lambert_rho = lambert_F * cot(pi/4 + lat/2)^lambert_n
    x = lambert_rho * sin(lambert_n * (long - ref_long))
    y = lambert_rho_ref - lambert_rho * cos(lambert_n * (long - ref_long))
    return (x, y)
end

# Determine partition shapes
println("Determining partition shapes...")
region_shapes = Dict{Int64, LibGEOS.MultiPolygon}()
partition_shapes = Dict{Int64, LibGEOS.GEOSGeom}()
for region_proto in region_list_proto.multi_polygon
    id = region_proto.id
    if !haskey(partition_list, id)
        continue
    end
    partition = partition_list[id]

    # Get region coordinates
    num_polygons = length(region_proto.polygon)
    region_coords = Vector{Vector{Vector{Vector{Float64}}}}()
    for (k, polygon_proto) in enumerate(region_proto.polygon)
        push!(region_coords, Vector{Vector{Vector{Float64}}}())
        for j in 1:(1 + length(polygon_proto.interior_border))
            if j == 1
                border_proto = polygon_proto.exterior_border
            else
                border_proto = polygon_proto.interior_border[j-1]
            end
            num_coords = length(border_proto.x)
            push!(region_coords[k], Vector{Vector{Float64}}())
            for i in 1:num_coords
                push!(region_coords[k][j], Vector{Float64}(2))
                region_coords[k][j][i][1] = border_proto.x[i]
                region_coords[k][j][i][2] = border_proto.y[i]
            end
        end
    end

    # Construct region shape
    region_shapes[id] = LibGEOS.MultiPolygon(region_coords)
    region_shape = region_shapes[id]

    # Add region shape to partition shape
    if haskey(partition_shapes, partition)
        partition_shapes[partition] = LibGEOS.union(partition_shapes[partition],
                                                    region_shape.ptr)
    else
        partition_shapes[partition] = LibGEOS.union(region_shape.ptr, region_shape.ptr)
    end
    
end

# Plot regions
println("Plotting...")
PyPlot.figure(figsize=(16, 16))
PyPlot.axis("off")
for (partition, shape) in partition_shapes
    partition_coords = GeoInterface.coordinates(LibGEOS.MultiPolygon(shape))
    for polygon_coords in partition_coords
        for border_coords in polygon_coords
            num_coords = length(border_coords)
            x = Vector{Float64}(num_coords)
            y = Vector{Float64}(num_coords)
            for i in 1:num_coords
                (x[i], y[i]) = lambert_projection(border_coords[i][1],
                                                  border_coords[i][2])
            end
            PyPlot.fill(x, y, color=color_list[partition])
            PyPlot.plot(x, y, "k-", linewidth=2.5)
        end
    end
end
for (region, shape) in region_shapes
    partition_coords = GeoInterface.coordinates(shape)
    for polygon_coords in partition_coords
        for border_coords in polygon_coords
            num_coords = length(border_coords)
            x = Vector{Float64}(num_coords)
            y = Vector{Float64}(num_coords)
            for i in 1:num_coords
                (x[i], y[i]) = lambert_projection(border_coords[i][1],
                                                  border_coords[i][2])
            end
            PyPlot.plot(x, y, "k-", linewidth=0.5)
        end
    end
end

# Export image
println("Exporting image...")
PyPlot.axis("tight")
PyPlot.axis("equal")
PyPlot.savefig(image_file)
