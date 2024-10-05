module Plot

import Base.Threads
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
            @inbounds border = polygon[1]  # External border
            @inbounds for j in 1:size(border, 2)
                @inbounds long::Float64 = border[1, j]
                @inbounds lat::Float64 = border[2, j]
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

struct Plotter
    county_boundaries::Dict{UInt, Vector{Vector{Array{Float64, 2}}}}
end

function Plotter(
    county_ids::Vector{UInt},
    )::Plotter

    # Load county boundaries from file
    county_boundaries = DataFiles.load_county_boundaries()

    # Filter unused counties
    county_ids_set = Set{UInt}(county_ids)
    county_boundaries = Dict(
        id => boundary
        for (id, boundary) in county_boundaries if id in county_ids_set)
    county_ids = [id for id in county_ids_set if haskey(county_boundaries, id)]
    sort!(county_ids)

    # Convert coordinates to Lambert projection
    lambert_projection = LambertProjection(county_boundaries)
    @Base.Threads.threads for county_id in county_ids
        multipolygon = county_boundaries[county_id]
        @inbounds for polygon in multipolygon
            @inbounds for border in polygon
                @inbounds for i in 1:size(border, 2)
                    @inbounds long = border[1,i]
                    @inbounds lat = border[2,i]
                    (x, y) = lambert_projection(long, lat)
                    @inbounds border[1,i] = x
                    @inbounds border[2,i] = y
                end
            end
        end
    end


    # Construct plotter
    return Plotter(county_boundaries)

end

function plot(plotter::Plotter)

    # Plot county borders
    @inbounds for multipolygon in values(plotter.county_boundaries)
        @inbounds for polygon in multipolygon
            @inbounds for border in polygon
                @inbounds x = border[1,:]
                @inbounds y = border[2,:]
                PyPlot.plot(x, y, "k-", linewidth=0.5)
            end
        end
    end

    # Show plot
    PyPlot.axis("off")
    PyPlot.axis("tight")
    PyPlot.axis("equal")
    PyPlot.show()

end

end  # module Plot
