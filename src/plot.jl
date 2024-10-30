module Plot

import Base.Threads
import ColorSchemes
import GLMakie
import Memoize

import ..Gerrymander
import ..DataFiles
import ..Geometry

@Memoize.memoize function wong_colorscheme()::Vector{GLMakie.RGBf}
    colors = [
        (86, 180, 233),
        (213, 94, 0),
        (0, 158, 115),
        (240, 228, 66),
        (0, 114, 178),
        (204, 121, 167),
        (230, 159, 0),
    ]
    return [GLMakie.RGBf(r/255, g/255, b/255) for (r, g, b) in colors]
end

function categorical_colors(data::Dict{UInt, UInt})::Dict{UInt, GLMakie.RGBf}
    colorscheme = wong_colorscheme()
    colors = Dict{UInt, GLMakie.RGBf}()
    sizehint!(colors, length(data))
    for (id, val) in data
        colors[id] = colorscheme[val % length(colorscheme) + 1]
    end
    return colors
end

function range_colors(
    data::Dict;
    log_scale::Bool = false,
    range_min::Union{Float64, Nothing} = nothing,
    range_max::Union{Float64, Nothing} = nothing,
    colorscheme::Union{ColorSchemes.ColorScheme, Nothing} = nothing,
    )::Dict{UInt, GLMakie.RGBf}

    # Find value range
    min_val::Float64 = Inf
    max_val::Float64 = -Inf
    for val::Float64 in values(data)
        min_val = min(min_val, val)
        max_val = max(max_val, val)
    end
    if isnothing(range_min)
        range_min = min_val
    end
    if isnothing(range_max)
        range_max = max_val
    end
    if log_scale
        range_min = log(range_min)
        range_max = log(range_max)
    end

    # Scale and shift to convert values to unit range
    shift = range_min
    scale = 1 / max(range_max - range_min, 1e-12)

    # Pick color from viridis colorscheme
    colors = Dict{UInt, GLMakie.RGBf}()
    sizehint!(colors, length(data))
    if isnothing(colorscheme)
        colorscheme = ColorSchemes.viridis
    end
    for (id, val::Float64) in data
        if log_scale
            val = log(val)
        end
        colors[id] = colorscheme[(val - shift) * scale]
    end

    return colors

end

mutable struct Plotter
    partitioner::Any
    county_ids::Vector{UInt}
    partition_ids::Vector{UInt}
    county_boundaries::Dict{UInt, Geometry.MultiPolygonCoords}
    boundary_lines::Vector{Geometry.PlotLine}
    frame_time::Float64  # sec

    function Plotter(
        partitioner::Any;
        frame_time::Float64 = 1.0,
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

        # Get county boundaries in Lambert projection
        county_boundaries = DataFiles.load_county_boundaries_lambert(county_ids)

        # Convert boundaries into lines
        boundary_lines = Geometry.county_boundaries_to_lines(county_boundaries)

        # Construct plotter
        return new(
            partitioner,
            county_ids,
            partition_ids,
            county_boundaries,
            boundary_lines,
            frame_time,
        )

    end

end

function plot(plotter::Plotter)

    # Compute partition shapes
    partition_shapes = Geometry.make_partition_shapes(
        plotter.county_boundaries,
        plotter.partitioner.partition_to_counties,
    )

    # Initialize plot
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], aspect=GLMakie.DataAspect())
    GLMakie.empty!(ax)
    GLMakie.hidespines!(ax)
    GLMakie.hidedecorations!(ax)

    # Plot partitions
    for partition_id in plotter.partition_ids
        color = pick_color(partition_id)
        @inbounds for polygon in partition_shapes[partition_id]
            GLMakie.poly!(
                polygon,
                color=GLMakie.RGBf(color[1], color[2], color[3]),
                stroke_depth_shift=0,
                strokecolor=:black,
                strokewidth=2,
            )
        end
    end

    # Plot county boundaries
    @inbounds for (x, y) in plotter.boundary_lines
        GLMakie.lines!(x, y, color = :black, linewidth = 0.25)
    end

    # Display plot
    GLMakie.display(fig)

    # Stall by waiting for user input
    println("Press enter to exit...")
    readline()

end

function animate!(plotter::Plotter)

    # Initialize plot
    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1, 1], aspect=GLMakie.DataAspect())
    GLMakie.empty!(ax)
    GLMakie.hidespines!(ax)
    GLMakie.hidedecorations!(ax)

    # State for frame updates
    step_time::Float64 = plotter.frame_time
    plot_time::Float64 = plotter.frame_time
    time_decay = 0.5
    is_paused::Bool = true
    display_is_paused::Bool = false
    show_counties::Bool = false

    function steps_per_frame()::Int
        frame_time = plotter.frame_time - (display_is_paused ? 0.0 : plot_time)
        num_steps = floor(Int, frame_time / step_time)
        num_steps = max(num_steps, 1)
    end

    function plot_partitions(colors::Dict{UInt, GLMakie.RGBf})

        # Compute partition shapes
        partition_shapes = Geometry.make_partition_shapes(
            plotter.county_boundaries,
            plotter.partitioner.partition_to_counties,
        )

        # Reset plot
        GLMakie.empty!(ax)

        # Plot partitions
        for partition_id in plotter.partition_ids
            color = colors[partition_id]
            for polygon in partition_shapes[partition_id]
                GLMakie.poly!(
                    polygon,
                    color=color,
                    stroke_depth_shift=0,
                    strokecolor=:black,
                    strokewidth=2,
                )
            end
        end

        # Plot county boundaries
        if show_counties
            @inbounds for (x, y) in plotter.boundary_lines
                GLMakie.lines!(x, y, color=:black, linewidth=0.25)
            end
        end

        # Display plot
        GLMakie.display(fig)

    end

    function plot_counties(colors::Dict{UInt, GLMakie.RGBf})

        # Compute partition shapes
        partition_shapes = Geometry.make_partition_shapes(
            plotter.county_boundaries,
            plotter.partitioner.partition_to_counties,
        )

        # Reset plot
        GLMakie.empty!(ax)

        # Plot counties
        for (county_id, multipolygon) in plotter.county_boundaries
            color = colors[county_id]
            for polygon in multipolygon
                plot_polygon = Vector{Vector{Tuple{Float64, Float64}}}()
                sizehint!(plot_polygon, length(polygon))
                for line in polygon
                    plot_line = Vector{Tuple{Float64, Float64}}()
                    sizehint!(plot_line, length(line))
                    push!(plot_polygon, plot_line)
                    for point in line
                        push!(plot_line, (point[1], point[2]))
                    end
                end
                GLMakie.poly!(
                    plot_polygon,
                    color=color,
                    stroke_depth_shift=0,
                    strokecolor=:black,
                    strokewidth=0.25,
                )
            end
        end

        # Plot partitions
        for partition_id in plotter.partition_ids
            @inbounds for polygon in partition_shapes[partition_id]
                GLMakie.poly!(
                    polygon,
                    color=:transparent,
                    stroke_depth_shift=0,
                    strokecolor=:black,
                    strokewidth=2,
                )
            end
        end

        # Display plot
        GLMakie.display(fig)

    end

    function plot_partition_data(
        data::Dict;
        log_scale::Bool = false,
        range_min::Union{Float64, Nothing} = nothing,
        range_max::Union{Float64, Nothing} = nothing,
        colorscheme::Union{ColorSchemes.ColorScheme, Nothing} = nothing,
        )
        colors = range_colors(
            data,
            log_scale=log_scale,
            range_min=range_min,
            range_max=range_max,
            colorscheme=colorscheme,
        )
        plot_partitions(colors)
    end

    function plot_county_data(
        data::Dict;
        log_scale::Bool = false,
        range_min::Union{Float64, Nothing} = nothing,
        range_max::Union{Float64, Nothing} = nothing,
        colorscheme::Union{ColorSchemes.ColorScheme, Nothing} = nothing,
        )
        colors = range_colors(
            data,
            log_scale=log_scale,
            range_min=range_min,
            range_max=range_max,
            colorscheme=colorscheme,
        )
        plot_counties(colors)
    end

    function update_plot()
        colors = categorical_colors(
            Dict{UInt, UInt}(id => id for id in plotter.partition_ids)
        )
        plot_partitions(colors)
    end

    # Initialize plot
    update_plot()

    function parse_command(command::String, params::String)
        if command == ""
        elseif command == "anim info" || command == "frame info"
            println()
            println("Animation info")
            println("--------------")
            println("Paused: ", is_paused)
            println("Paused display: ", display_is_paused)
            println("Show counties: ", show_counties)
            println("Partitioner step time: ", step_time * 1e3, " ms")
            println("Plot time: ", plot_time * 1e3, " ms")
            println("Frame time: ", plotter.frame_time * 1e3, " ms")
            println("Steps per frame: ", steps_per_frame())
            println()
        elseif command == "pause" || command == "unpause"
            if isempty(params)
                is_paused = !is_paused
            else
                is_paused = parse(Bool, params)
            end
            if is_paused
                println("Pausing...")
            else
                println("Unpausing...")
            end
        elseif command == "pause plot"
            if isempty(params)
                display_is_paused = !display_is_paused
            else
                display_is_paused = parse(Bool, params)
            end
            if display_is_paused
                println("Pausing display...")
            else
                println("Unpausing display...")
            end
        elseif command == "reset plot" || command == "reset anim" || command == "update plot"
            update_plot()
        elseif command == "show counties"
            show_counties = parse(Bool, params)
            update_plot()
        elseif command == "save image" || command == "save plot"
            println("Saving image to $params")
            GLMakie.save(params, fig, px_per_unit=8)
        elseif hasfield(typeof(plotter.partitioner), :parse_command_func)
            plotter.partitioner.parse_command_func(
                command,
                params,
                plot_partition_data,
                plot_county_data,
            )
        else
            println("Unrecognized command: ", command)
        end
    end

    # Spawn thread to monitor commands from user
    loop_is_active::Base.Threads.Atomic{Bool} = Base.Threads.Atomic{Bool}(true)
    command_lock = ReentrantLock()
    command_channel = Channel{String}() do ch
        lock(command_lock)
        while loop_is_active[]
            unlock(command_lock)
            print("Command: ")
            commands = readline()
            if strip(commands) != ""
                put!(ch, commands)
            end
            lock(command_lock)
        end
        unlock(command_lock)
    end

    # Animation loop
    while loop_is_active[]
        loop_start_time = time()

        if !is_paused

            # Perform partitioner steps
            num_steps = steps_per_frame()
            for _ in 1:num_steps
                Gerrymander.step!(plotter.partitioner)
            end
            step_end_time = time()
            step_time *= (1 - time_decay)
            step_time += time_decay * (step_end_time - loop_start_time) / num_steps

            # Update plot
            if !display_is_paused
                update_plot()
                plot_end_time = time()
                plot_time *= (1 - time_decay)
                plot_time += time_decay * (plot_end_time - step_end_time)
            end

        end

        # Process user commands
        if !isopen(command_channel)
            loop_is_active[] = false
        elseif isready(command_channel)
            lock(command_lock)
            commands = take!(command_channel)
            for command in split(commands, ";")
                command = lowercase(command)
                command_split = split(command, "=", limit=2)
                command = String(strip(command_split[1]))
                params = length(command_split) > 1 ? String(strip(command_split[2])) : ""
                if command == "exit"
                    println("Exiting...")
                    loop_is_active[] = false
                    break
                else
                    try
                        parse_command(command, params)
                    catch err
                        println("Error while processing command \"$command\"")
                        println(err)
                    end
                end
            end
            unlock(command_lock)
        end
        if !loop_is_active[]
            break
        end

        # Sleep to achieve desired frame rate
        sleep_start_time = time()
        sleep_time = plotter.frame_time - (sleep_start_time - loop_start_time)
        if sleep_time > 0
            sleep(sleep_time)
        end

    end

    # Clean up
    GLMakie.closeall()

end

end  # module Plot
