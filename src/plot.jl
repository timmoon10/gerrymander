module Plot

import Base.Threads
import DataStructures
import Memoize
import PyPlot
import ..DataFiles

"Parameters for Lambert conformal conic projection"
struct LambertProjection
    ref_long::Float64
    ref_lat::Float64
    n::Float64
    F::Float64
    rho_ref::Float64
end

function LambertProjection(
    county_boundaries::Dict{UInt, Vector{Vector{Array{Float64, 2}}}},
    )::LambertProjection

    # Determine region boundaries
    min_long::Float64 = 0.0
    max_long::Float64 = -180.0
    min_lat::Float64 = 90.0
    max_lat::Float64 = -90.0
    @inbounds for multipolygon in values(county_boundaries)
        @inbounds for polygon in multipolygon
            @inbounds line = polygon[1]  # External border
            @inbounds for j in 1:size(line, 2)
                @inbounds long::Float64 = line[1, j]
                @inbounds lat::Float64 = line[2, j]
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
    county_boundaries::Dict{UInt, Vector{Vector{Array{Float64, 2}}}},
    )::Vector{Tuple{Vector{Float64}, Vector{Float64}}}

    "Convert boundaries into line segments"
    function make_segments(
        county_boundaries::Dict{UInt, Vector{Vector{Array{Float64, 2}}}},
        )::Dict{Tuple{Float64, Float64}, Set{Tuple{Float64, Float64}}}
        segments = Dict{Tuple{Float64, Float64}, Set{Tuple{Float64, Float64}}}()
        @inbounds for multipolygon in values(county_boundaries)
            @inbounds for polygon in multipolygon
                @inbounds for line in polygon

                    # Skip if current line is empty
                    num_segments = size(line, 2) - 1
                    if num_segments < 1
                        continue
                    end

                    # Break down line into line segments
                    @inbounds coord1 = (line[1, 1], line[2, 1])
                    if !haskey(segments, coord1)
                        segments[coord1] = Set{Tuple{Float64, Float64}}()
                    end
                    @inbounds for i in 1:num_segments
                        @inbounds coord2 = (line[1, i+1], line[2, i+1])
                        if !haskey(segments, coord2)
                            segments[coord2] = Set{Tuple{Float64, Float64}}()
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
        segments::Dict{Tuple{Float64, Float64}, Set{Tuple{Float64, Float64}}},
        coord::Tuple{Float64, Float64},
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

struct Plotter
    county_boundaries::Dict{UInt, Vector{Vector{Array{Float64, 2}}}}
    boundary_lines::Vector{Tuple{Vector{Float64}, Vector{Float64}}}
end

function Plotter(
    partitions::Dict{UInt, UInt},
    )::Plotter

    # Get lists of counties and partitions
    county_ids = Set{UInt}()
    partition_ids = Set{UInt}()
    for (county_id, partition_id) in partitions
        push!(county_ids, county_id)
        push!(partition_ids, partition_id)
    end
    county_ids = collect(county_ids)
    partition_ids = collect(partition_ids)
    sort!(county_ids)
    sort!(partition_ids)

    # Get county boundaries
    county_boundaries = Dict{UInt, Vector{Vector{Array{Float64, 2}}}}()
    let
        full_county_boundaries = DataFiles.load_county_boundaries()
        for county_id in county_ids
            county_boundaries[county_id] = full_county_boundaries[county_id]
        end
    end

    # Convert coordinates to Lambert projection
    lambert_projection = LambertProjection(county_boundaries)
    @Base.Threads.threads for county_id in county_ids
        multipolygon = county_boundaries[county_id]
        @inbounds for polygon in multipolygon
            @inbounds for line in polygon
                @inbounds for i in 1:size(line, 2)
                    @inbounds long = line[1,i]
                    @inbounds lat = line[2,i]
                    (x, y) = lambert_projection(long, lat)
                    @inbounds line[1,i] = x
                    @inbounds line[2,i] = y
                end
            end
        end
    end

    # Convert boundaries into lines
    boundary_lines = county_boundaries_to_lines(county_boundaries)

    # Construct plotter
    return Plotter(county_boundaries, boundary_lines)

end

function plot(plotter::Plotter)

    # Plot county boundaries
    @inbounds for (x, y) in plotter.boundary_lines
        PyPlot.plot(x, y, "k-", linewidth=0.5)
    end

    # Show plot
    PyPlot.axis("off")
    PyPlot.axis("tight")
    PyPlot.axis("equal")
    PyPlot.show()

end

end  # module Plot
