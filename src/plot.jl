module Plot

import Base.GC
import Base.Threads
import DataStructures
import GeoInterface
import LibGEOS
import Memoize
import PyCall

using PyCall
@pyimport matplotlib.animation as animation
import PyPlot

import ..DataFiles
using ..Constants: MultiPolygonCoords, PolygonCoords
import ..SimulatedAnnealing  ### TODO Remove

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
    return [(r/255, g/255, b/255) for (r, g, b) in colors]
end

function pick_color(i::UInt)::Tuple{Float64, Float64, Float64}
    colors = color_list()
    return colors[i % length(colors) + 1]
end

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
    county_ids = collect(keys(county_boundaries))
    coord_bounds = Vector{Tuple{Float64, Float64, Float64, Float64}}(
        undef,
        length(county_ids),
    )
    @Base.Threads.threads for i in 1:length(county_ids)
        county_id = county_ids[i]
        coord_bounds[i] = multipolygon_bounds(county_boundaries[county_id])
    end
    (min_long, max_long, min_lat, max_lat) = coord_bounds[1]
    @inbounds for bounds in coord_bounds[2:end]
        @inbounds min_long = min(min_long, bounds[1])
        @inbounds max_long = max(max_long, bounds[2])
        @inbounds min_lat = min(min_lat, bounds[3])
        @inbounds max_lat = max(max_lat, bounds[4])
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

# x- and y-coordinates for a line
PlotLine = Tuple{Vector{Float64}, Vector{Float64}}

mutable struct Plotter
    partitioner::Any
    county_ids::Vector{UInt}
    partition_ids::Vector{UInt}
    county_boundaries::Dict{UInt, MultiPolygonCoords}
    boundary_lines::Vector{PlotLine}
end

function Plotter(
    partitioner::Any,
    )::Plotter

    # Get lists of counties and partitions
    county_ids = Set{UInt}()
    partition_ids = Set{UInt}()
    for (county_id, partition_id) in partitioner.county_to_partition
        push!(county_ids, county_id)
        push!(partition_ids, partition_id)
    end
    county_ids = collect(county_ids)
    partition_ids = collect(partition_ids)
    sort!(county_ids)
    sort!(partition_ids)

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

    # Compute grid size for downsampling
    coord_bounds = Vector{Tuple{Float64, Float64, Float64, Float64}}(
        undef,
        length(county_ids),
    )
    @Base.Threads.threads for i in 1:length(county_ids)
        county_id = county_ids[i]
        coord_bounds[i] = multipolygon_bounds(county_boundaries[county_id])
    end
    (min_x, max_x, min_y, max_y) = coord_bounds[1]
    @inbounds for bounds in coord_bounds[2:end]
        @inbounds min_x = min(min_x, bounds[1])
        @inbounds max_x = max(max_x, bounds[2])
        @inbounds min_y = min(min_y, bounds[3])
        @inbounds max_y = max(max_y, bounds[4])
    end
    grid_size = min(max_x - min_x, max_y - min_y) / 1024

    # Downsample boundaries
    downsample_county_boundaries!(county_boundaries, grid_size)

    # Convert boundaries into lines
    boundary_lines = county_boundaries_to_lines(county_boundaries)

    # Construct plotter
    out = Plotter(
        partitioner,
        county_ids,
        partition_ids,
        county_boundaries,
        boundary_lines,
    )
    return out

end

function make_partition_shapes(plotter)::Dict{UInt, Vector{Vector{PlotLine}}}

    # Objects from plotter
    partition_ids = plotter.partition_ids
    partition_to_counties = plotter.partitioner.partition_to_counties
    county_boundaries = plotter.county_boundaries

    # Iterate through partitions
    partition_shapes = Vector{Vector{Vector{PlotLine}}}(undef, length(partition_ids))
    @Base.Threads.threads for i in 1:length(partition_ids)

        # Get county shapes
        multipolygon = MultiPolygonCoords()
        @inbounds for county_id in partition_to_counties[partition_ids[i]]
            append!(multipolygon, county_boundaries[county_id])
        end

        # Compute union of county shapes
        multipolygon = LibGEOS.MultiPolygon(multipolygon)
        multipolygon = LibGEOS.unaryUnion(multipolygon)
        multipolygon = LibGEOS.MultiPolygon(multipolygon)
        multipolygon = GeoInterface.coordinates(multipolygon)

        # Convert partition shape to lines
        plot_multipolygon = Vector{Vector{PlotLine}}()
        sizehint!(plot_multipolygon, length(multipolygon))
        @inbounds for polygon in multipolygon
            plot_polygon = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
            push!(plot_multipolygon, plot_polygon)
            sizehint!(plot_polygon, length(polygon))
            @inbounds for line in polygon
                x = Vector{Float64}()
                y = Vector{Float64}()
                push!(plot_polygon, (x, y))
                sizehint!(x, length(line))
                sizehint!(y, length(line))
                @inbounds for point in line
                    @inbounds push!(x, point[1])
                    @inbounds push!(y, point[2])
                end
            end
        end
        partition_shapes[i] = plot_multipolygon

    end
    partition_shapes = Dict{UInt, Vector{Vector{PlotLine}}}(
        partition_ids[i] => shape for (i, shape) in enumerate(partition_shapes))

    return partition_shapes

end

function plot(plotter::Plotter)

    # Compute partition shapes
    partition_shapes = make_partition_shapes(plotter)

    # Plot partitions
    for partition_id in plotter.partition_ids
        color = pick_color(partition_id)
        @inbounds for polygon in partition_shapes[partition_id]
            @inbounds for (x, y) in polygon
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

function animate!(
    plotter::Plotter;
    frame_interval::UInt = UInt(1000),
    )

    # Initialize plot
    fig = PyPlot.figure()
    ax = PyPlot.axes(xlim=(-0.05, 0.05), ylim=(-0.05, 0.05))
    PyPlot.axis("off")
    PyPlot.axis("tight")
    PyPlot.axis("equal")

    function plot_frame()

        # Compute partition shapes
        partition_shapes = make_partition_shapes(plotter)

        # Reset plot
        ax.clear()

        # Plot partitions
        for partition_id in plotter.partition_ids
            color = pick_color(partition_id)
            @inbounds for polygon in partition_shapes[partition_id]
                @inbounds for (x, y) in polygon
                    ax.fill(x, y, color=color)
                    ax.plot(x, y, "k-", linewidth=1)
                end
            end
        end

        # Plot county boundaries
        @inbounds for (x, y) in plotter.boundary_lines
            ax.plot(x, y, "k-", linewidth=0.25)
        end

        # Plot options
        PyPlot.axis("off")
        PyPlot.axis("tight")
        PyPlot.axis("equal")

    end

    # State for frame updates
    lock = ReentrantLock()
    steps_per_frame::Int = 1
    step_time::Float64 = frame_interval
    plot_time::Float64 = frame_interval
    running_mean_decay = 0.5

    function update_frame(frame::Int)

        # Skip frame if earlier frame is late
        if !trylock(lock)
            return
        end

        # Work around strange segfaults from PyCall
        Base.GC.enable(false)

        # Perform partitioner steps
        step_start_time = time()
        for i in 1:steps_per_frame
            SimulatedAnnealing.step!(plotter.partitioner)
        end
        step_end_time = time()

        # Plot partitioner state
        plot_frame()
        plot_end_time = time()

        # Update running averages for run times
        step_time *= (1 - running_mean_decay)
        step_time += (
            running_mean_decay * 1e3 * (step_end_time - step_start_time)
            / steps_per_frame
        )
        plot_time *= (1 - running_mean_decay)
        plot_time += (
            running_mean_decay * 1e3 * (plot_end_time - step_end_time)
        )

        # Update number of steps per frame
        steps_per_frame = floor(
            Int,
            (frame_interval * 0.8 - plot_time) / step_time,
        )
        steps_per_frame = max(steps_per_frame, 1)

        # Clean up
        Base.GC.enable(true)
        unlock(lock)

    end

    # Start animation
    anim = animation.FuncAnimation(
        fig,
        update_frame,
        init_func=plot_frame,
        interval=frame_interval,
        cache_frame_data=false,
    )
    PyPlot.show()

end

end  # module Plot
