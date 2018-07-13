#!/usr/bin/julia
#
# Plot county partition
#
using ProtoBuf
using PyCall
import DataStructures
using PyPlot
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Import partition data
println("Reading partition data...")
part_list = Dict{UInt32, UInt32}()
num_parts = 0
part_data = readdlm(partition_file, '\t')
for row in 1:size(part_data, 1)
    geoid = UInt32(part_data[row, 1])
    part = UInt32(part_data[row, 2])
    part_list[geoid] = part
    num_parts = max(part, num_parts)
end

# Generate color scheme
println("Generating plot colors...")
push!(PyVector(pyimport("sys")["path"]), project_dir * "/randomcolor-py")
@pyimport randomcolor
color_list = []
for i = 1:num_parts
    color = randomcolor.RandomColor()[:generate]()[1]
    push!(color_list, color)
end

# Read boundary data from protobuf
println("Reading boundary data...")
bounds_list_proto = CountyBoundariesList()
open(county_bounds_file, "r") do f
    readproto(f, bounds_list_proto)
end

# Initialize Lambert conformal conic projection
println("Determining projection parameters...")
min_long = 0
max_long = -180
min_lat = 90
max_lat = -90
for bounds_proto in bounds_list_proto.county_bounds
    geoid = bounds_proto.geoid
    if haskey(part_list, geoid)
        for polygon_proto in bounds_proto.polygon
            for long in polygon_proto.x
                min_long = min(long, min_long)
                max_long = max(long, max_long)
            end
            for lat in polygon_proto.y
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
    long = long_deg * pi / 180
    lat  = lat_deg * pi / 180
    lambert_rho = lambert_F * cot(pi/4 + lat/2)^lambert_n
    x = lambert_rho * sin(lambert_n * (long - ref_long))
    y = lambert_rho_ref - lambert_rho * cos(lambert_n * (long - ref_long))
    return (x, y)
end

# Determine internal and external segments
println("Determining partition borders...")
exterior_segments = Dict{Tuple{Float32, Float32, Float32, Float32}, UInt32}()
interior_segments = Set{Tuple{Float32, Float32, Float32, Float32}}()
for bounds_proto in bounds_list_proto.county_bounds
    geoid = bounds_proto.geoid
    if !haskey(part_list, geoid)
        continue
    end
    part = part_list[geoid]
    for polygon_proto in bounds_proto.polygon
        num_coords = length(polygon_proto.x) - 1
        for i in 1:num_coords
            segment = (polygon_proto.x[i], polygon_proto.y[i],
                       polygon_proto.x[i+1], polygon_proto.y[i+1])
            reverse_segment = (segment[3], segment[4],
                               segment[1], segment[2])
            if haskey(exterior_segments, segment)
                if exterior_segments[segment] == part
                    delete!(exterior_segments, segment)
                    push!(interior_segments, reverse_segment)
                end
            else
                exterior_segments[reverse_segment] = part
            end
        end
    end
end

# Plot counties
println("Plotting counties...")
PyPlot.figure(figsize=(16, 16))
PyPlot.axis("off")
for bounds_proto in bounds_list_proto.county_bounds
    geoid = bounds_proto.geoid
    if !haskey(part_list, geoid)
        continue
    end
    part = part_list[geoid]
    color = color_list[part]
    for polygon_proto in bounds_proto.polygon
        num_points = length(polygon_proto.x)
        x = Array{Float32, 1}(num_points)
        y = Array{Float32, 1}(num_points)
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
