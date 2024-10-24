module Plot

import GLMakie
import Memoize

import ..Gerrymander
import ..DataFiles
import ..Geometry

@Memoize.memoize function color_list()::Vector{GLMakie.RGBf}
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

function pick_color(i::UInt)::GLMakie.RGBf
    colors = color_list()
    return colors[i % length(colors) + 1]
end

mutable struct Plotter
    partitioner::Any
    county_ids::Vector{UInt}
    partition_ids::Vector{UInt}
    county_boundaries::Dict{UInt, Geometry.MultiPolygonCoords}
    boundary_lines::Vector{Geometry.PlotLine}
    frame_time::Float64  # sec
    is_paused::Bool

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
            0.125,
            false,
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

    function steps_per_frame()::Int
        num_steps = floor(Int, (plotter.frame_time - plot_time) / step_time)
        num_steps = max(num_steps, 1)
    end

    function update_plot()

        # Compute partition shapes
        partition_shapes = Geometry.make_partition_shapes(
            plotter.county_boundaries,
            plotter.partitioner.partition_to_counties,
        )

        # Reset plot
        GLMakie.empty!(ax)

        # Plot partitions
        for partition_id in plotter.partition_ids
            color = pick_color(partition_id)
            @inbounds for polygon in partition_shapes[partition_id]
                @inbounds for (x, y) in polygon
                    GLMakie.poly!(
                        polygon,
                        color=color,
                        stroke_depth_shift=0,
                        strokecolor=:black,
                        strokewidth=2,
                    )
                end
            end
        end

        # Plot county boundaries
        @inbounds for (x, y) in plotter.boundary_lines
            GLMakie.lines!(x, y, color=:black, linewidth=0.25)
        end

        # Display plot
        GLMakie.display(fig)

    end

    # Initialize plot
    update_plot()

    function parse_command(command::String)
        if command == "pause" || command == "unpause"
            if plotter.is_paused
                println("Unpausing...")
            else
                println("Pausing...")
            end
            plotter.is_paused = !plotter.is_paused
        elseif command == "anim info" || command == "frame info"
            println()
            println("Animation info")
            println("--------------")
            println("Paused: ", plotter.is_paused)
            println("Partitioner step time: ", step_time * 1e3, " ms")
            println("Plot time: ", plot_time * 1e3, " ms")
            println("Frame time: ", plotter.frame_time * 1e3, " ms")
            println("Steps per frame: ", steps_per_frame())
            println()
        elseif command == "reset plot" || command == "reset anim"
            update_plot()
        elseif hasfield(typeof(plotter.partitioner), :parse_command_func)
            plotter.partitioner.parse_command_func(command)
            update_plot()
        else
            println("Unrecognized command: ", command)
        end
    end

    # Spawn thread to monitor commands from user
    command_lock = ReentrantLock()
    command_channel = Channel{String}() do ch
        loop_is_active = true
        while loop_is_active
            lock(command_lock)
            unlock(command_lock)
            print("Command: ")
            for command in split(readline(), ";")
                command = lowercase(strip(command))
                if command == "exit" || command == "quit"
                    put!(ch, "exit")
                    loop_is_active = false
                    break
                end
                if command != ""
                    put!(ch, command)
                end
            end
        end
    end

    # Animation loop
    loop_is_active = true
    while loop_is_active
        loop_start_time = time()

        if !plotter.is_paused

            # Perform partitioner steps
            num_steps = steps_per_frame()
            for _ in 1:num_steps
                Gerrymander.step!(plotter.partitioner)
            end
            step_end_time = time()

            # Update plot
            update_plot()
            plot_end_time = time()

            # Update times
            step_time *= (1 - time_decay)
            step_time += time_decay * (step_end_time - loop_start_time) / num_steps
            plot_time *= (1 - time_decay)
            plot_time += time_decay * (plot_end_time - step_end_time)

        end

        # Process user commands
        if !isopen(command_channel)
            loop_is_active = false
        elseif isready(command_channel)
            lock(command_lock)
            command = take!(command_channel)
            if command == "exit"
                println("Exiting...")
                loop_is_active = false
            else
                parse_command(command)
            end
            unlock(command_lock)
        end
        if !loop_is_active
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
