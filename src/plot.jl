module Plot

import Base.Threads
import DataStructures
import GeoInterface
import LibGEOS
import Memoize
import PyPlot

import ..DataFiles
using ..Constants: MultiPolygonCoords

@Memoize.memoize function color_list()::Vector{Tuple{Float64, Float64, Float64}}
    colors = [
        (86, 180, 233),
        (213, 94, 0),
        (0, 158, 115),
        (240, 228, 66),
        (0, 114, 178),
        (204, 121, 167),
        (230, 159, 0),
    ]
    return [(r/256, g/256, b/256) for (r, g, b) in colors]
end

function hash_color(x)::Tuple{Float64, Float64, Float64}
    colors = color_list()
    return colors[hash(x) % length(colors) + 1]
end

"Parameters for Lambert conformal conic projection"
struct LambertProjection
    ref_long::Float64
    ref_lat::Float64
    n::Float64
    F::Float64
    rho_ref::Float64
end

function LambertProjection(
    county_boundaries::Dict{UInt, MultiPolygonCoords},
    )::LambertProjection

    # Determine region boundaries
    min_long::Float64 = 0.0
    max_long::Float64 = -180.0
    min_lat::Float64 = 90.0
    max_lat::Float64 = -90.0
    @inbounds for multipolygon in values(county_boundaries)
        @inbounds for polygon in multipolygon
            @inbounds for coord in polygon[1]
                @inbounds long::Float64 = coord[1]
                @inbounds lat::Float64 = coord[2]
                min_long = min(long, min_long)
                max_long = max(long, max_long)
                min_lat = min(lat, min_lat)
                max_lat = max(lat, max_lat)
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

function (lambert_projection::LambertProjection)(
    long_deg::Float64,
    lat_deg::Float64,
    )::Tuple{Float64, Float64}
    ref_long = lambert_projection.ref_long
    ref_lat = lambert_projection.ref_lat
    n = lambert_projection.n
    F = lambert_projection.F
    rho_ref = lambert_projection.rho_ref
    deg_to_rad::Float64 = pi / 180
    long = deg_to_rad * long_deg
    lat = deg_to_rad * lat_deg
    rho = F * cot(pi/4 + lat/2)^n
    x = rho * sin(n * (long - ref_long))
    y = rho_ref - rho * cos(n * (long - ref_long))
    return (x, y)
end

function county_boundaries_to_lines(
    county_boundaries::Dict{UInt, MultiPolygonCoords},
    )::Vector{Tuple{Vector{Float64}, Vector{Float64}}}

    "Convert boundaries into line segments"
    function make_segments(
        county_boundaries::Dict{UInt, MultiPolygonCoords},
        )::Dict{Vector{Float64}, Set{Vector{Float64}}}
        segments = Dict{Vector{Float64}, Set{Vector{Float64}}}()
        @inbounds for multipolygon in values(county_boundaries)
            @inbounds for polygon in multipolygon
                @inbounds for line in polygon

                    # Skip if current line is empty
                    num_segments = length(line) - 1
                    if num_segments < 1
                        continue
                    end

                    # Break down line into line segments
                    @inbounds coord1 = line[1]
                    if !haskey(segments, coord1)
                        segments[coord1] = Set{Vector{Float64}}()
                    end
                    @inbounds for i in 1:num_segments
                        @inbounds coord2 = line[i+1]
                        if !haskey(segments, coord2)
                            segments[coord2] = Set{Vector{Vector{Float64}}}()
                        end
                        push!(segments[coord1], coord2)
                        push!(segments[coord2], coord1)
                        coord1 = coord2
                    end

                end
            end
        end
        return segments
    end

    "Traverse line segments, removing as we go"
    function walk_segments!(
        segments::Dict{Vector{Float64}, Set{Vector{Float64}}},
        coord::Vector{Float64},
        )::Tuple{Vector{Float64}, Vector{Float64}}

        # Walk over segments until stuck
        xs = [coord[1]]
        ys = [coord[2]]
        while haskey(segments, coord)

            # Walk over segment and delete
            neighbors = segments[coord]
            neighbor = pop!(neighbors)
            if isempty(neighbors)
                delete!(segments, coord)
            end
            neighbor_neighbors = segments[neighbor]
            delete!(neighbor_neighbors, coord)
            if isempty(neighbor_neighbors)
                delete!(segments, neighbor)
            end
            coord = neighbor

            # Record new coordinate
            push!(xs, coord[1])
            push!(ys, coord[2])

        end
        return (xs, ys)

    end

    # Convert county boundaries into lines
    lines = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    segments = make_segments(county_boundaries)
    while !isempty(segments)
        coord = first(keys(segments))
        (xs1, ys1) = walk_segments!(segments, coord)
        (xs2, ys2) = walk_segments!(segments, coord)
        reverse!(xs2)
        reverse!(ys2)
        push!(lines, ([xs2; xs1[2:end]], [ys2; ys1[2:end]]))
    end
    return lines

end

mutable struct Plotter
    county_ids::Vector{UInt}
    partition_ids::Vector{UInt}
    partition_to_counties::Dict{UInt, Set{UInt}}
    county_boundaries::Dict{UInt, MultiPolygonCoords}
    boundary_lines::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
    county_shapes::Dict{UInt, LibGEOS.MultiPolygon}
    partition_shapes::Dict{UInt, LibGEOS.MultiPolygon}
end

function Plotter(
    county_to_partition::Dict{UInt, UInt},
    )::Plotter

    # Get lists of counties and partitions
    county_ids = Set{UInt}()
    partition_ids = Set{UInt}()
    for (county_id, partition_id) in county_to_partition
        push!(county_ids, county_id)
        push!(partition_ids, partition_id)
    end
    county_ids = collect(county_ids)
    partition_ids = collect(partition_ids)
    sort!(county_ids)
    sort!(partition_ids)

    # Counties in each partition
    partition_to_counties = Dict{UInt, Set{UInt}}(
        partition_id => Set{UInt}() for partition_id in partition_ids)
    for (county_id, partition_id) in county_to_partition
        push!(partition_to_counties[partition_id], county_id)
    end

    # Get county boundaries
    county_boundaries = Dict{UInt, MultiPolygonCoords}()
    let
        full_county_boundaries = DataFiles.load_county_boundaries()
        @inbounds for county_id in county_ids
            county_boundaries[county_id] = full_county_boundaries[county_id]
        end
    end

    "Apply Lambert projection in-place to multipolygon coordinates"
    function multipolygon_to_lambert!(
        multipolygon::MultiPolygonCoords,
        lambert_projection::LambertProjection,
        )
        @inbounds for polygon in multipolygon
            @inbounds for line in polygon
                @inbounds for coord in line
                    @inbounds long = coord[1]
                    @inbounds lat = coord[2]
                    (x, y) = lambert_projection(long, lat)
                    @inbounds coord[1] = x
                    @inbounds coord[2] = y
                end
            end
        end
    end

    # Convert coordinates to Lambert projection
    lambert_projection = LambertProjection(county_boundaries)
    @Base.Threads.threads for county_id in county_ids
        multipolygon_to_lambert!(
            county_boundaries[county_id],
            lambert_projection,
        )
    end

    # Convert boundaries into lines
    boundary_lines = county_boundaries_to_lines(county_boundaries)

    # Construct plotter
    out = Plotter(
        county_ids,
        partition_ids,
        partition_to_counties,
        county_boundaries,
        boundary_lines,
        Dict{UInt, LibGEOS.MultiPolygon}(),
        Dict{UInt, LibGEOS.MultiPolygon}(),
    )
    reset_shapes!(out, reset_counties=true, reset_partitions=true)
    return out

end

function reset_shapes!(
    plotter::Plotter;
    reset_counties::Bool=true,
    reset_partitions::Bool=true,
    )::Plotter

    # Objects from plotter
    county_ids = plotter.county_ids
    partition_ids = plotter.partition_ids
    partition_to_counties = plotter.partition_to_counties
    county_boundaries = plotter.county_boundaries

    # County shapes
    if reset_counties
        county_shapes = Vector{LibGEOS.MultiPolygon}(undef, length(county_ids))
        @Base.Threads.threads for i in 1:length(county_ids)
            multipolygon = county_boundaries[county_ids[i]]
            county_shapes[i] = LibGEOS.MultiPolygon(multipolygon)
        end
        county_shapes = Dict{UInt, LibGEOS.MultiPolygon}(
            county_ids[i] => shape for (i, shape) in enumerate(county_shapes))
        plotter.county_shapes = county_shapes
    end

    # Partition shapes
    if reset_partitions
        partition_shapes = Vector{LibGEOS.MultiPolygon}(undef, length(partition_ids))
        @Base.Threads.threads for i in 1:length(partition_ids)
            multipolygon = MultiPolygonCoords()
            @inbounds for county_id in partition_to_counties[partition_ids[i]]
                append!(multipolygon, county_boundaries[county_id])
            end
            multipolygon = LibGEOS.MultiPolygon(multipolygon)
            multipolygon = LibGEOS.unaryUnion(multipolygon)
            multipolygon = LibGEOS.MultiPolygon(multipolygon)
            partition_shapes[i] = multipolygon
        end
        partition_shapes = Dict{UInt, LibGEOS.MultiPolygon}(
            partition_ids[i] => shape for (i, shape) in enumerate(partition_shapes))
        plotter.partition_shapes = partition_shapes
    end

    return plotter

end

function plot(plotter::Plotter)

    # Plot partitions
    for (partition_id, shape) in plotter.partition_shapes
        color = hash_color(partition_id)
        multipolygon = GeoInterface.coordinates(shape)
        @inbounds for polygon in multipolygon
            @inbounds for line in polygon
                x = Vector{Float64}(undef, length(line))
                y = Vector{Float64}(undef, length(line))
                @inbounds for (point_id, point) in enumerate(line)
                    @inbounds x[point_id] = point[1]
                    @inbounds y[point_id] = point[2]
                end
                PyPlot.fill(x, y, color=color)
                PyPlot.plot(x, y, "k-", linewidth=1)
            end
        end
    end

    # Plot county boundaries
    @inbounds for (x, y) in plotter.boundary_lines
        PyPlot.plot(x, y, "k-", linewidth=0.25)
    end

    # Show plot
    PyPlot.axis("off")
    PyPlot.axis("tight")
    PyPlot.axis("equal")
    PyPlot.show()

end

end  # module Plot
