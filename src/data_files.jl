module DataFiles

import DataStructures
import DelimitedFiles
import Downloads
import JSON
import Memoize

import ..Constants

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
    file_base = "tl_2020_" * lpad(state_id, 2, "0") * "_tabblock20"
    download_dir = joinpath(root_dir(), "data", "downloads")
    json_file = joinpath(download_dir, file_base * ".json")
    if !isfile(json_file)
        println("Downloading census block data for " * state_names()[state_id] * "...")
        census_url = "https://www2.census.gov/geo/tiger/TIGER2020/TABBLOCK20/"
        zip_file = joinpath(download_dir, file_base * ".zip")
        shape_file = joinpath(download_dir, file_base * ".shp")
        Downloads.download(census_url * file_base * ".zip", zip_file)
        run(`unzip $zip_file -d $download_dir`)
        run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
    end
    return json_file
end

function maybe_download_geography_data()::String
    file_base = "tl_2020_us_county"
    download_dir = joinpath(root_dir(), "data", "downloads")
    json_file = joinpath(download_dir, file_base * ".json")
    if !isfile(json_file)
        println("Downloading geography data...")
        census_url = "https://www2.census.gov/geo/tiger/TIGER2020/COUNTY/"
        zip_file = joinpath(download_dir, file_base * ".zip")
        shape_file = joinpath(download_dir, file_base * ".shp")
        Downloads.download(census_url * file_base * ".zip", zip_file)
        run(`unzip $zip_file -d $download_dir`)
        run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
    end
    return json_file
end

function maybe_parse_county_populations(state_ids::AbstractVector{UInt})::String

    # Return immediately if county population data file exists
    county_populations_file = joinpath(
        root_dir(),
        "results",
        "county_populations.tsv",
    )
    if isfile(county_populations_file)
        return county_populations_file
    end

    # Get census block data
    json_files = [
        maybe_download_census_block_data(state_id)
        for state_id in state_ids]

    function parse_state_data(
        json_file::String,
        state_id::UInt,
        )::Array{Any, 2}

        # Constants
        deg_to_rad::Float64 = pi / 180

        # Parse JSON file to get position and population of census blocks
        sums = Dict{UInt, Array{Float64}}()
        for block_json::Dict{String, Any} in JSON.parsefile(json_file)["features"]

            # Approximate census block position
            sum_x::Float64 = 0.0
            sum_y::Float64 = 0.0
            sum_z::Float64 = 0.0
            num_coords::UInt = 0
            geometry = block_json["geometry"]["type"]
            if geometry == "Polygon"
                for polygon::Vector{Any} in block_json["geometry"]["coordinates"]
                    for coord::Vector{Float64} in polygon
                        long::Float64 = deg_to_rad * coord[1]
                        lat::Float64 = deg_to_rad * coord[2]
                        sum_x += cos(long) * sin(lat)
                        sum_y += sin(long) * sin(lat)
                        sum_z += cos(lat)
                        num_coords += 1
                    end
                end
            elseif geometry == "MultiPolygon"
                for multipolygon::Vector{Any} in block_json["geometry"]["coordinates"]
                    for polygon::Vector{Any} in multipolygon
                        for coord::Vector{Float64} in polygon
                            long::Float64 = deg_to_rad * coord[1]
                            lat::Float64 = deg_to_rad * coord[2]
                            sum_x += cos(long) * sin(lat)
                            sum_y += sin(long) * sin(lat)
                            sum_z += cos(lat)
                            num_coords += 1
                        end
                    end
                end
            end
            x::Float64 = Constants.earth_radius * sum_x / num_coords
            y::Float64 = Constants.earth_radius * sum_y / num_coords
            z::Float64 = Constants.earth_radius * sum_z / num_coords

            # Initialize county data
            props_json = block_json["properties"]
            county_id = parse(UInt, props_json["COUNTYFP20"])
            geoid::UInt = state_id * 1000 + county_id
            if !haskey(sums, geoid)
                sums[geoid] = zeros(Float64, 7)
            end
            county_sums = sums[geoid]
            pop::Float64 = props_json["POP20"]

            # Record population and position of census block
            county_sums[1] += pop * x
            county_sums[2] += pop * y
            county_sums[3] += pop * z
            county_sums[4] += pop * x * x
            county_sums[5] += pop * y * y
            county_sums[6] += pop * z * z
            county_sums[7] += pop

        end

        # Compute population, votes, and position of counties
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
            state_county_data[i, 2] = pop
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
    println("Exporting county data...")
    DelimitedFiles.writedlm(county_populations_file, county_data, '\t')
    return county_populations_file

end


end  # module DataFiles
