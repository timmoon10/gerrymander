#!/usr/bin/julia
#
# Plot county partition
#
using ProtoBuf
import DataStructures
using PyPlot
include(AbstractString(dirname(@__FILE__)) * "/common.jl")

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

# Generate color scheme with https://github.com/kevinwuhoo/randomcolor-py
println("Generating plot colors...")
color_list = []
for i = 1:num_parts
    color_list = [color_list; chomp(readall(`$color_generator`))]
end

# Read boundary data from protobuf
println("Reading boundary data...")
bounds_list_proto = CountyBoundariesList()
open(county_bounds_file, "r") do f
    readproto(f, bounds_list_proto)
end

# Determine projection parameters
# Note: Lambert conformal conic projection
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
ref_long = (min_long + max_long) / 2
ref_lat = (min_lat + max_lat) / 2
lambert_n = (log(cos(max_lat) * sec(min_lat))
             / log(tan(pi/4 + min_lat/2) * cot(pi/4 + max_lat/2)))
lambert_F = (cos(max_lat) * tan(pi/4+max_lat/2)^lambert_n)/ lambert_n
lambert_rho_ref = lambert_F * cot(pi/4 + ref_lat/2)^lambert_n

# Plotting each county
println("Plotting counties...")
PyPlot.figure(figsize=(16, 16))
PyPlot.axis("off")
for bounds_proto in bounds_list_proto.county_bounds

    # Get county properties
    geoid = bounds_proto.geoid
    if !haskey(part_list, geoid)
        continue
    end
    part = part_list[geoid]
    color = color_list[part]

    # Plot each polygon in county
    for polygon_proto in bounds_proto.polygon
        
        # Lambert conformal conic projection
        num_points = length(polygon_proto.x)
        x = Array{Float32, 1}(num_points)
        y = Array{Float32, 1}(num_points)
        for i = 1:num_points
            long = polygon_proto.x[i] * pi / 180
            lat = polygon_proto.y[i] * pi / 180
            lambert_rho = lambert_F * cot(pi/4 + lat/2)^lambert_n
            x[i] = lambert_rho * sin(lambert_n * (long - ref_long))
            y[i] = lambert_rho_ref - lambert_rho * cos(lambert_n * (long - ref_long))
        end

        # Plot polygon
        PyPlot.fill(x, y, color=color)
        PyPlot.plot(x, y, "k-", linewidth=0.5)
        
    end
    
end
PyPlot.axis("tight")
PyPlot.axis("equal")
PyPlot.savefig(image_file)
