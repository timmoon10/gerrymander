module DataFiles

import Base.Threads
import DataStructures
import DelimitedFiles
import Downloads
import GDAL_jll
import JSON
import Memoize
import Serialization

import ..Constants
import ..Geometry
using ..Geometry: MultiPolygonCoords, PolygonCoords

@Memoize.memoize function root_dir()::String
    return realpath(dirname(realpath(@__DIR__)))
end

@Memoize.memoize function state_names()::DataStructures.OrderedDict{UInt, String}
    (data, _) = DelimitedFiles.readdlm(
        joinpath(root_dir(), "data", "state_codes.tsv"),
        '\t',
        header=true,
        comments=true,
    )
    names = DataStructures.OrderedDict{UInt, String}()
    for row in 1:size(data, 1)
        names[data[row, 3]] = data[row, 1]
    end
    return names
end

function maybe_download_census_block_data(state_id::UInt)::String

    # Return immediately if census block data is ready
    file_base = "tl_2020_" * lpad(state_id, 2, "0") * "_tabblock20"
    download_dir = joinpath(root_dir(), "data", "downloads")
    json_file = joinpath(download_dir, file_base * ".json")
    if isfile(json_file)
        return json_file
    end

    # Download data
    println("Downloading census block data for " * state_names()[state_id] * "...")
    census_url = "https://www2.census.gov/geo/tiger/TIGER2020/TABBLOCK20/"
    zip_file = joinpath(download_dir, file_base * ".zip")
    Downloads.download(census_url * file_base * ".zip", zip_file)

    # Extract and parse shape file
    shape_file = joinpath(download_dir, file_base * ".shp")
    run(`unzip $zip_file -d $download_dir`)
    run(`$(GDAL_jll.ogr2ogr_path()) -f "GeoJSON" $json_file $shape_file`)

    return json_file

end

function maybe_download_geography_data()::String

    # Return immediately if census block data is ready
    file_base = "tl_2020_us_county"
    download_dir = joinpath(root_dir(), "data", "downloads")
    json_file = joinpath(download_dir, file_base * ".json")
    if isfile(json_file)
        return json_file
    end

    # Download data
    println("Downloading geography data...")
    census_url = "https://www2.census.gov/geo/tiger/TIGER2020/COUNTY/"
    zip_file = joinpath(download_dir, file_base * ".zip")
    Downloads.download(census_url * file_base * ".zip", zip_file)

    # Extract and parse shape file
    shape_file = joinpath(download_dir, file_base * ".shp")
    run(`unzip $zip_file -d $download_dir`)
    run(`$(GDAL_jll.ogr2ogr_path()) -f "GeoJSON" $json_file $shape_file`)

    return json_file
end

"Load county populations from file

Each row corresponds to a county and has 6 fields: county ID,
population, x-coordinate, y-coordinate, z-coordinate, variance.

"
function load_county_populations(state_ids::AbstractVector{UInt})::Array{Any, 2}

    # Return immediately if cached data exists
    cache_file = joinpath(root_dir(), "results", "county_populations.tsv")
    if isfile(cache_file)
        (county_data, _) = DelimitedFiles.readdlm(
            cache_file,
            '\t',
            header=true,
        )
        return county_data
    end

    # Get census block data
    json_files = [maybe_download_census_block_data(state_id) for state_id in state_ids]

    "Compute sum of coordinates on unit sphere"
    function coords_sum(
        coords::Vector{Any},
        )::Tuple{Float64, Float64, Float64}
        deg_to_rad::Float64 = pi / 180
        sx::Float64 = 0
        sy::Float64 = 0
        sz::Float64 = 0
        @inbounds for coord::Vector{Float64} in coords
            @inbounds long::Float64 = deg_to_rad * coord[1]
            @inbounds lat::Float64 = deg_to_rad * coord[2]
            @inbounds (sin_long, cos_long) = sincos(long)
            @inbounds (sin_lat, cos_lat) = sincos(lat)
            sx += cos_long * sin_lat
            sy += sin_long * sin_lat
            sz += cos_lat
        end
        return (sx, sy, sz)
    end

    "Parse population data for a state"
    function parse_state_data(
        json_file::String,
        state_id::UInt,
        )::Array{Any, 2}

        # Parse JSON file
        state_data = JSON.parsefile(json_file)["features"]

        # Parse census block data in parallel
        make_sums = () -> zeros(Float64, 7)
        thread_sums = [
            DataStructures.DefaultDict{UInt, Array{Float64}}(make_sums)
            for _ in 1:Base.Threads.nthreads()]
        @Base.Threads.threads for i in 1:length(state_data)
            block_data = state_data[i]

            # Approximate census block position
            sum_x::Float64 = 0.0
            sum_y::Float64 = 0.0
            sum_z::Float64 = 0.0
            num_coords::UInt = 0
            geometry = block_data["geometry"]["type"]
            if geometry == "Polygon"
                for polygon::Vector{Any} in block_data["geometry"]["coordinates"]
                    (sx, sy, sz) = coords_sum(polygon)
                    sum_x += sx
                    sum_y += sy
                    sum_z += sz
                    num_coords += length(polygon)
                end
            elseif geometry == "MultiPolygon"
                for multipolygon::Vector{Any} in block_data["geometry"]["coordinates"]
                    for polygon::Vector{Any} in multipolygon
                        (sx, sy, sz) = coords_sum(polygon)
                        sum_x += sx
                        sum_y += sy
                        sum_z += sz
                        num_coords += length(polygon)
                    end
                end
            end
            x::Float64 = Constants.earth_radius * sum_x / num_coords
            y::Float64 = Constants.earth_radius * sum_y / num_coords
            z::Float64 = Constants.earth_radius * sum_z / num_coords

            # Census block properties
            block_props = block_data["properties"]
            pop::Float64 = block_props["POP20"]
            county_id = parse(UInt, block_props["COUNTYFP20"])
            geoid::UInt = state_id * 1000 + county_id

            # Accumulate sums
            county_sums = thread_sums[Base.Threads.threadid()][geoid]
            county_sums[1] += pop * x
            county_sums[2] += pop * y
            county_sums[3] += pop * z
            county_sums[4] += pop * x * x
            county_sums[5] += pop * y * y
            county_sums[6] += pop * z * z
            county_sums[7] += pop

        end

        # Reduce sums across parallel threads
        sums = DataStructures.DefaultDict{UInt, Array{Float64}}(make_sums)
        for thread_sums_i in thread_sums
            for (geoid, county_sums) in thread_sums_i
                sums[geoid] += county_sums
            end
        end

        # Compute population and position of counties
        state_county_data = Array{Any, 2}(undef, length(sums), 6)
        for (i, (geoid, county_sums)) in enumerate(sums)
            pop = county_sums[7]
            x = county_sums[1] / pop
            y = county_sums[2] / pop
            z = county_sums[3] / pop
            var_x = county_sums[4] / pop - x * x
            var_y = county_sums[5] / pop - y * y
            var_z = county_sums[6] / pop - z * z
            var = max(var_x + var_y + var_z, 0.0)
            state_county_data[i, 1] = geoid
            state_county_data[i, 2] = UInt(pop)
            state_county_data[i, 3] = x
            state_county_data[i, 4] = y
            state_county_data[i, 5] = z
            state_county_data[i, 6] = var
        end
        return state_county_data

    end

    # Compute county population data for each state
    println("Computing county populations...")
    county_data = ["geoid" "pop" "pos_x" "pos_y" "pos_z" "pos_var"]
    for (json_file, state_id) in zip(json_files, state_ids)
        county_data = [county_data; parse_state_data(json_file, state_id)]
    end

    # Output results to file
    DelimitedFiles.writedlm(cache_file, county_data, '\t')
    println("Saved county populations at $cache_file...")
    return county_data[2:end, :]

end

function county_boundaries_to_file(
    county_boundaries::Dict{UInt, MultiPolygonCoords},
    file::String,
    )::Nothing

    # Convert county boundary data to represent lines as 2D arrays
    # Note: Representing lines as vector of vectors results in many
    # dynamic dispatches.
    county_ids = collect(keys(county_boundaries))
    file_county_boundaries = Vector{Vector{Vector{Array{Float64, 2}}}}(
        undef,
        length(county_ids),
    )
    @Base.Threads.threads for i in 1:length(county_ids)
        multipolygon = county_boundaries[county_ids[i]]
        file_multipolygon = Vector{Vector{Array{Float64, 2}}}()
        sizehint!(file_multipolygon, length(multipolygon))
        @inbounds file_county_boundaries[i] = file_multipolygon
        @inbounds for (polygon_id, polygon) in enumerate(multipolygon)
            file_polygon = Vector{Array{Float64, 2}}()
            sizehint!(file_polygon, length(polygon))
            push!(file_multipolygon, file_polygon)
            @inbounds for (line_id, line) in enumerate(polygon)
                file_line = Array{Float64, 2}(undef, (2, length(line)))
                push!(file_polygon, file_line)
                @inbounds for (point_id, point) in enumerate(line)
                    @inbounds file_line[:, point_id] = point
                end
            end
        end
    end
    file_county_boundaries = Dict{UInt, Vector{Vector{Array{Float64, 2}}}}(
        id => file_county_boundaries[i] for (i, id) in enumerate(county_ids))

    # Save data to file
    Serialization.serialize(file, file_county_boundaries)

end

function county_boundaries_from_file(
    file::String,
    )::Dict{UInt, MultiPolygonCoords}

    # Load data from file
    file_county_boundaries = Serialization.deserialize(file)

    # Convert county boundary data to expected format
    county_ids = collect(keys(file_county_boundaries))
    county_boundaries = Vector{MultiPolygonCoords}(undef, length(county_ids))
    @Base.Threads.threads for i in 1:length(county_ids)
        file_multipolygon = file_county_boundaries[county_ids[i]]
        multipolygon = MultiPolygonCoords()
        sizehint!(multipolygon, length(file_multipolygon))
        @inbounds county_boundaries[i] = multipolygon
        @inbounds for (polygon_id, file_polygon) in enumerate(file_multipolygon)
            polygon = PolygonCoords()
            sizehint!(polygon, length(file_polygon))
            push!(multipolygon, polygon)
            @inbounds for (line_id, file_line) in enumerate(file_polygon)
                line = Vector{Vector{Float64}}()
                sizehint!(line, size(file_line, 2))
                push!(polygon, line)
                @inbounds for point_id in 1:size(file_line, 2)
                    @inbounds push!(line, file_line[:, point_id])
                end
            end
        end
    end
    county_boundaries = Dict{UInt, MultiPolygonCoords}(
        id => county_boundaries[i] for (i, id) in enumerate(county_ids))

    return county_boundaries

end

function load_county_boundaries()::Dict{UInt, MultiPolygonCoords}

    # Return immediately if county boundary data file exists
    cache_file = joinpath(root_dir(), "results", "county_boundaries.bin")
    if isfile(cache_file)
        return county_boundaries_from_file(cache_file)
    end

    # Get geography data
    json_file = maybe_download_geography_data()
    geography_data = JSON.parsefile(json_file)["features"]

    "Parse polygon boundary"
    function parse_polygon(polygon_data::Vector{Any})::PolygonCoords
        polygon = Vector{Vector{Vector{Float64}}}(undef, length(polygon_data))
        @inbounds for (line_id, line_data::Vector{Any}) in enumerate(polygon_data)
            @inbounds polygon[line_id] = [coord for coord in line_data]
        end
        return polygon
    end

    # Parse county boundaries
    println("Parsing county boundaries...")
    county_boundaries = Dict{UInt, MultiPolygonCoords}()
    for county_data in geography_data
        id = parse(UInt, county_data["properties"]["GEOID"])
        geometry_type = county_data["geometry"]["type"]
        if geometry_type == "Polygon"
            polygon_data = county_data["geometry"]["coordinates"]
            county_boundaries[id] = [parse_polygon(polygon_data)]
        elseif geometry_type == "MultiPolygon"
            multipolygon_data = county_data["geometry"]["coordinates"]
            county_boundaries[id] = [
                parse_polygon(polygon_data)
                for polygon_data in multipolygon_data]
        end
    end

    # Write results to file
    county_boundaries_to_file(county_boundaries, cache_file)
    println("Saved county boundaries at $cache_file...")

    return county_boundaries

end

function load_county_boundaries_lambert(
    county_ids;
    resolution::Int = 1024,
    )::Dict{UInt, MultiPolygonCoords}

    # Return immediately if cache file exists
    cache_file = joinpath(root_dir(), "results", "county_boundaries_lambert.bin")
    if isfile(cache_file)
        return county_boundaries_from_file(cache_file)
    end

    # List of counties
    county_ids::Vector{UInt} = collect(county_ids)

    # Get county boundaries
    county_boundaries = Dict{UInt, MultiPolygonCoords}()
    let
        full_county_boundaries = load_county_boundaries()
        @inbounds for county_id in county_ids
            county_boundaries[county_id] = full_county_boundaries[county_id]
        end
    end

    # Convert coordinates to Lambert projection
    lambert_projection = Geometry.LambertProjection(county_boundaries)
    Geometry.apply_lambert!(county_boundaries, lambert_projection)

    # Downsample boundaries if needed
    if resolution > 0

        # Compute grid size for downsampling
        coord_bounds = Vector{Tuple{Float64, Float64, Float64, Float64}}(
            undef,
            length(county_ids),
        )
        @Base.Threads.threads for i in 1:length(county_ids)
            county_id = county_ids[i]
            coord_bounds[i] = Geometry.multipolygon_bounds(county_boundaries[county_id])
        end
        (min_x, max_x, min_y, max_y) = coord_bounds[1]
        @inbounds for (min_xi, max_xi, min_yi, max_yi) in coord_bounds[2:end]
            min_x = min(min_x, min_xi)
            max_x = max(max_x, max_xi)
            min_y = min(min_y, min_yi)
            max_y = max(max_y, max_yi)
        end
        grid_size = min(max_x - min_x, max_y - min_y) / resolution

        # Downsample boundaries
        Geometry.downsample_county_boundaries!(county_boundaries, grid_size)

    end

    # Write results to file
    county_boundaries_to_file(county_boundaries, cache_file)
    println("Saved county boundaries (in Lambert projection) at $cache_file...")

    return county_boundaries

end

end  # module DataFiles
