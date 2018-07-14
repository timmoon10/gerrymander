#!/usr/bin/julia
#
# Plot partitions
#
using PyCall
import ProtoBuf, LibGEOS, GeoInterface, PyPlot
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Import partition data
println("Importing partition data...")
part_data = readdlm(partition_file, '\t')
num_regions = size(part_data, 1)
part_list = Dict{Int64, Int64}()
num_parts = 0
for row in 1:num_regions
    geoid = Int64(part_data[row, 1])
    part = Int64(part_data[row, 2])
    part_list[geoid] = part
    num_parts = max(part, num_parts)
end

# Generate color scheme
println("Generating plot colors...")
push!(PyCall.PyVector(pyimport("sys")["path"]), project_dir * "/randomcolor-py")
@pyimport randomcolor
color_list = []
for i = 1:num_parts
    color = randomcolor.RandomColor()[:generate]()[1]
    push!(color_list, color)
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
    if haskey(part_list, id)
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
    if !haskey(part_list, id)
        continue
    end
    partition = part_list[id]

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
println("Plotting regions...")
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
            PyPlot.plot(x, y, "k-", linewidth=4.0)
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

exit()

# Plot regions
println("Plotting regions...")
PyPlot.figure(figsize=(16, 16))
PyPlot.axis("off")
region_shapes = Dict{Int64, LibGEOS.MultiPolygon}()
partition_shapes = Dict{Int64, LibGEOS.MultiPolygon}()


for bounds_proto in bounds_list_proto.county_bounds
    geoid = bounds_proto.geoid
    if !haskey(part_list, geoid)
        continue
    end
    part = part_list[geoid]
    color = color_list[part]
    for polygon_proto in bounds_proto.polygon
        num_points = length(polygon_proto.x)
        x = Vector{Float64}(num_points)
        y = Vector{Float64}(num_points)
        for i = 1:num_points
            (x[i], y[i]) = lambert_projection(polygon_proto.x[i],
                                              polygon_proto.y[i])
        end
        PyPlot.fill(x, y, color=color)
        PyPlot.plot(x, y, "k-", linewidth=0.5)
    end
end
for (long1, lat1, long2, lat2) in keys(exterior_segments)
    (x1, y1) = lambert_projection(long1, lat1)
    (x2, y2) = lambert_projection(long2, lat2)
    PyPlot.plot([x1, x2], [y1, y2], "k-", linewidth=4.0)
end

# Export image
println("Exporting image...")
PyPlot.axis("tight")
PyPlot.axis("equal")
PyPlot.savefig(image_file)
