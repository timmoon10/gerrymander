module Plot

import Base.GC
import GLMakie
import Memoize

import ..Gerrymander
import ..DataFiles
import ..Geometry

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

mutable struct Plotter
    partitioner::Any
    county_ids::Vector{UInt}
    partition_ids::Vector{UInt}
    county_boundaries::Dict{UInt, Geometry.MultiPolygonCoords}
    boundary_lines::Vector{Geometry.PlotLine}
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

function animate!(
    plotter::Plotter;
    frame_interval::UInt = UInt(1000),
    )::Nothing

    # Initialize plot
    fig = PyPlot.figure()
    ax = PyPlot.axes()
    PyPlot.axis("off")

    function plot_frame()

        # Compute partition shapes
        partition_shapes = Geometry.make_partition_shapes(
            plotter.county_boundaries,
            plotter.partitioner.partition_to_counties,
        )

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
    anim_lock = ReentrantLock()
    steps_per_frame::Int = 1
    step_time::Float64 = frame_interval
    plot_time::Float64 = frame_interval
    running_mean_decay = 0.5

    "Update animation frame"
    function update_frame(frame::Int)::Nothing

        # Skip frame if paused
        if plotter.is_paused
            return
        end

        # Skip frame if mutex is locked
        if !trylock(anim_lock)
            return
        end

        # Work around strange segfaults from PyCall
        Base.GC.enable(false)

        # Perform partitioner steps
        step_start_time = time()
        for i in 1:steps_per_frame
            Gerrymander.step!(plotter.partitioner)
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
        unlock(anim_lock)

    end

    # Start animation
    anim = animation.FuncAnimation(
        fig,
        update_frame,
        init_func=plot_frame,
        interval=frame_interval,
        cache_frame_data=false,
    )

    "Logic for key presses"
    function on_key(event)::Nothing

        # Lock mutex to avoid interfering with animation
        lock(anim_lock)

        # Work around strange segfaults from PyCall
        Base.GC.enable(false)

        # Basic key press logic
        if event.key == "escape"
            exit()
        end
        if event.key == "p"
            plotter.is_paused = !plotter.is_paused
        end

        # Partitioner key press logic
        if hasfield(typeof(plotter.partitioner), :on_key_func)
            plotter.partitioner.on_key_func(event)
        end

        # Clean up
        Base.GC.enable(true)
        unlock(anim_lock)

    end

    # Register logic for key presses
    fig.canvas.mpl_connect("key_press_event", on_key)

    # Start animation
    PyPlot.show()

end

end  # module Plot
