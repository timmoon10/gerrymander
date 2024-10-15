module Geometry

import Base.Threads
import GeoInterface
import LibGEOS

using ..Constants: MultiPolygonCoords, PolygonCoords

function multipolygon_bounds(
    multipolygon::MultiPolygonCoords,
    )::Tuple{Float64, Float64, Float64, Float64}
    min_x::Float64 = Inf
    max_x::Float64 = -Inf
    min_y::Float64 = Inf
    max_y::Float64 = -Inf
    @inbounds for polygon in multipolygon
        @inbounds for point in polygon[1]
            @inbounds x = point[1]
            @inbounds y = point[2]
            min_x = min(x, min_x)
            max_x = max(x, max_x)
            min_y = min(y, min_y)
            max_y = max(y, max_y)
        end
    end
    return (min_x, max_x, min_y, max_y)
end

function downsample_county_boundaries!(
    county_boundaries::Dict{UInt, MultiPolygonCoords},
    grid_size::Float64,
    )

    # Helper functions to convert to/from grid
    to_grid = let grid_size::Float64 = grid_size
        function to_grid(x::Float64)::Int
            return div(x, grid_size)
        end
        to_grid
    end
    from_grid = let grid_size::Float64 = grid_size
        function from_grid(i::Int)::Float64
            return i * grid_size
        end
        from_grid
    end

    # Find grid points in each county
    county_ids = collect(keys(county_boundaries))
    county_grid_bounds = Dict{UInt, Tuple{Int, Int, Int, Int}}(
        county_id => (0, 0, 0, 0) for county_id in county_ids)
    county_grid_points = Dict{UInt, Set{Tuple{Int, Int}}}(
        county_id => Set{Tuple{Int, Int}}() for county_id in county_ids)
    @Base.Threads.threads for county_id in county_ids

        # Grid bounds
        multipolygon = county_boundaries[county_id]
        (min_x, max_x, min_y, max_y) = multipolygon_bounds(multipolygon)
        min_i = 2 * div(to_grid(min_x), 2) - 2
        max_i = 2 * div(to_grid(max_x), 2) + 4
        min_j = 2 * div(to_grid(min_y), 2) - 2
        max_j = 2 * div(to_grid(max_y), 2) + 4
        county_grid_bounds[county_id]  = (min_i, max_i, min_j, max_j)

        # Find grid points in multipolygon
        # Note: Sample at half resolution
        multipolygon = LibGEOS.MultiPolygon(multipolygon)
        grid_points = county_grid_points[county_id]
        sizehint!(grid_points, (max_i-min_i+1) * (max_j-min_j+1))
        for i in min_i:2:max_i
            x = from_grid(i)
            for j in min_j:2:max_j
                y = from_grid(j)
                point = LibGEOS.Point(x, y)
                if LibGEOS.contains(multipolygon, point)
                    push!(grid_points, (i, j))
                end
            end
        end

    end

    # Find county for each grid point
    grid_point_counties = Dict{Tuple{Int, Int}, Int}()
    num_points = sum(
        length(grid_points)
        for grid_points in values(county_grid_points))
    sizehint!(grid_point_counties, num_points)
    for (county_id, grid_points) in county_grid_points
        for grid_point in grid_points
            grid_point_counties[grid_point] = county_id
        end
    end

    # Helper function to make full grid square
    function make_full_grid(i::Int, j::Int)::PolygonCoords
        x1 = from_grid(i)
        x2 = from_grid(i+1)
        y1 = from_grid(j)
        y2 = from_grid(j+1)
        return [[[x1, y1], [x1, y2], [x2, y2], [x2, y1], [x1, y1]]]
    end

    # Helper function to make half grid triangle
    function make_half_grid(orientation::Int, i::Int, j::Int)::PolygonCoords
        x1 = from_grid(i)
        x2 = from_grid(i+1)
        y1 = from_grid(j)
        y2 = from_grid(j+1)
        if orientation == 1
            return [[[x1, y1], [x1, y2], [x2, y1], [x1, y1]]]
        elseif orientation == 2
            return [[[x1, y2], [x2, y2], [x1, y1], [x1, y2]]]
        elseif orientation == 3
            return [[[x2, y2], [x2, y1], [x1, y2], [x2, y2]]]
        elseif orientation == 4
            return [[[x2, y1], [x1, y1], [x2, y2], [x2, y1]]]
        else
            return []
        end
    end

    # Downsample each county boundary
    @Base.Threads.threads for county_id in county_ids

        # Iterate over half-resolution grid
        multipolygon = MultiPolygonCoords()
        (min_i, max_i, min_j, max_j) = county_grid_bounds[county_id]
        for i in min_i:2:max_i-1
            for j in min_j:2:max_j-1

                # Check if corners are in multipolygon
                corner_counties = (
                    get(grid_point_counties, (i,j), -1),
                    get(grid_point_counties, (i,j+2), -1),
                    get(grid_point_counties, (i+2,j+2), -1),
                    get(grid_point_counties, (i+2,j), -1),
                )
                corner_counties = (
                    corner_counties[1],
                    corner_counties[2],
                    corner_counties[3],
                    corner_counties[4],
                    corner_counties[1],
                    corner_counties[2],
                    corner_counties[3],
                )

                # Process grid squares in groups of 4
                corner_is = (i, i, i+1, i+1)
                corner_js = (j, j+1, j+1, j)
                for corner in 1:4
                    corner1_id = corner_counties[corner]
                    corner2_id = corner_counties[corner+1]
                    corner3_id = corner_counties[corner+2]
                    corner4_id = corner_counties[corner+3]
                    corner_i = corner_is[corner]
                    corner_j = corner_js[corner]
                    if corner1_id == county_id
                        if (
                            corner2_id != county_id
                            && corner2_id == corner3_id
                            && corner3_id == corner4_id
                            && corner4_id == corner2_id
                            )
                            push!(
                                multipolygon,
                                make_half_grid(corner, corner_i, corner_j),
                            )
                        else
                            push!(
                                multipolygon,
                                make_full_grid(corner_i, corner_j),
                            )
                        end
                    elseif (
                        corner2_id == county_id
                        && corner3_id == county_id
                        && corner4_id == county_id
                        )
                        orientation = corner <= 2 ? corner + 2 : corner - 2
                        push!(
                            multipolygon,
                            make_half_grid(orientation, corner_i, corner_j),
                        )
                    end
                end

            end
        end

        # Compute union of downsampled multipolygon
        if !isempty(multipolygon)
            multipolygon = LibGEOS.MultiPolygon(multipolygon)
            multipolygon = LibGEOS.unaryUnion(multipolygon)
            multipolygon = LibGEOS.MultiPolygon(multipolygon)
            multipolygon = GeoInterface.coordinates(multipolygon)
        end

        # Update county boundaries with downsampled multipolygon
        county_boundaries[county_id] = multipolygon

    end

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

end  # module Geometry
