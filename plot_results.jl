#
# Plot partitions
#
import DelimitedFiles
import GeoInterface
import LibGEOS
import PyPlot
import Serialization
include(joinpath(dirname(@__FILE__), "common.jl"))

# Import partition data
println("Importing partition data...")
partition_data = DelimitedFiles.readdlm(partition_file, '\t')
partitions = Dict{Int64, Int64}()
for row in 1:size(partition_data, 1)
    geoid = partition_data[row, 1]
    partition = partition_data[row, 2]
    partitions[geoid] = partition
end
max_partition = maximum(values(partitions))

# Generate color scheme
println("Generating plot colors...")
color_list = [
    (86,180,233), (213,94,0), (0,158,115), (240,228,66),
    (0,114,178), (204,121,167), (230,159,0),
]
color_list = [(r/256, g/256, b/256) for (r,g,b) in color_list]

# Read boundary data from file
println("Importing geography data...")
regions = Serialization.deserialize(geography_data_file)

# Struct with parameters for Lambert projection
struct LambertProjection
    ref_long::Float64
    ref_lat::Float64
    n::Float64
    F::Float64
    rho_ref::Float64
end

# Constructor to initialize parameters for Lambert projection
function LambertProjection(
    regions::Dict{Int64, Vector{Vector{Array{Float64, 2}}}},
    partitions::Dict{Int64, Int64},
    )::LambertProjection

    # Determine region boundaries
    min_long::Float64 = 0.0
    max_long::Float64 = -180.0
    min_lat::Float64 = 90.0
    max_lat::Float64 = -90.0
    for (id, region) in regions
        if haskey(partitions, id)
            for polygon in region
                border = polygon[1] # External border
                for j in 1:size(border, 2)
                    long::Float64 = border[1, j]
                    lat::Float64 = border[2, j]
                    min_long = min(long, min_long)
                    max_long = max(long, max_long)
                    min_lat = min(lat, min_lat)
                    max_lat = max(lat, max_lat)
                end
            end
        end
    end
    min_long *= pi / 180
    max_long *= pi / 180
    min_lat *= pi / 180
    max_lat *= pi / 180

    # Compute Lambert projection parameters
    ref_long = (min_long + max_long) / 2
    ref_lat = (min_lat + max_lat) / 2
    n = (
        log(cos(max_lat) * sec(min_lat))
        / log(tan(pi/4 + min_lat/2) * cot(pi/4 + max_lat/2))
    )
    F = (cos(max_lat) * tan(pi/4+max_lat/2)^n) / n
    rho_ref = F * cot(pi/4 + ref_lat/2)^n
    return LambertProjection(ref_long, ref_lat, n, F, rho_ref)

end

# Functions to apply Lambert projection
function (lambert_projection::LambertProjection)(
    long_deg::Float64,
    lat_deg::Float64,
    )::Tuple{Float64, Float64}
    ref_long = lambert_projection.ref_long
    ref_lat = lambert_projection.ref_lat
    n = lambert_projection.n
    F = lambert_projection.F
    rho_ref = lambert_projection.rho_ref
    long = pi / 180 * long_deg
    lat = pi / 180 * lat_deg
    rho = F * cot(pi/4 + lat/2)^n
    x = rho * sin(n * (long - ref_long))
    y = rho_ref - rho * cos(n * (long - ref_long))
    return (x, y)
end
function (lambert_projection::LambertProjection)(
    coords::Vector{Vector{Float64}},
    )::Array{Float64, 2}
    num_coords = length(coords)
    proj_coords = Array{Float64, 2}(undef, (2, num_coords))
    for i in 1:num_coords
        coord = coords[i]
        (x, y) = lambert_projection(coord[1], coord[2])
        proj_coords[1,i] = x
        proj_coords[2,i] = y
    end
    return proj_coords
end

# Initialize Lambert conformal conic projection
println("Initializing map projection...")
lambert_projection = LambertProjection(regions, partitions)

# Function to construct LibGEOS shapes for each region and partition
function construct_shapes(
    regions::Dict{Int64, Vector{Vector{Array{Float64, 2}}}},
    partitions::Dict{Int64, Int64},
    )::Tuple{Dict{Int64, LibGEOS.MultiPolygon}, Dict{Int64, LibGEOS.GEOSGeom}}
    region_shapes = Dict{Int64, LibGEOS.MultiPolygon}()
    partition_shapes = Dict{Int64, LibGEOS.GEOSGeom}()

    # Construct LibGEOS shape for each region
    for (id, region) in regions
        if !haskey(partitions, id)
            continue
        end
        partition::Int64 = partitions[id]

        # Get region coordinates
        region_coords = Vector{Vector{Vector{Vector{Float64}}}}()
        sizehint!(region_coords, length(region))
        for polygon in region
            push!(region_coords, Vector{Vector{Vector{Float64}}}())
            sizehint!(region_coords[end], length(polygon))
            for border in polygon
                push!(region_coords[end], Vector{Vector{Float64}}())
                sizehint!(region_coords[end][end], size(border, 2))
                for j in 1:size(border, 2)
                    push!(region_coords[end][end], border[:, j])
                end
            end
        end

        # Construct region shape
        region_shapes[id] = LibGEOS.MultiPolygon(region_coords)
        region_shape = region_shapes[id]

        # Add region shape to partition shape
        if haskey(partition_shapes, partition)
            partition_shapes[partition] = LibGEOS.union(
                partition_shapes[partition],
                region_shape.ptr,
            )
        else
            partition_shapes[partition] = LibGEOS.union(
                region_shape.ptr,
                region_shape.ptr,
            )
        end

    end

    return (region_shapes, partition_shapes)
end

# Construct shapes
println("Constructing partition shapes...")
(region_shapes, partition_shapes) = construct_shapes(regions, partitions)

# Plot regions
println("Plotting...")
PyPlot.figure(figsize=(16, 16))
PyPlot.axis("off")
for (partition, shape) in partition_shapes
    partition_coords = GeoInterface.coordinates(LibGEOS.MultiPolygon(shape))
    for polygon_coords in partition_coords
        for border_coords in polygon_coords
            plot_coords = lambert_projection(border_coords)
            x = plot_coords[1,:]
            y = plot_coords[2,:]
            color_id = (partition % length(color_list)) + 1
            PyPlot.fill(x, y, color=color_list[color_id])
            PyPlot.plot(x, y, "k-", linewidth=2.5)
        end
    end
end
for shape in values(region_shapes)
    region_coords = GeoInterface.coordinates(shape)
    for polygon_coords in region_coords
        for border_coords in polygon_coords
            plot_coords = lambert_projection(border_coords)
            x = plot_coords[1,:]
            y = plot_coords[2,:]
            PyPlot.plot(x, y, "k-", linewidth=0.5)
        end
    end
end

# Export image
println("Exporting image...")
PyPlot.axis("tight")
PyPlot.axis("equal")
PyPlot.savefig(image_file)
