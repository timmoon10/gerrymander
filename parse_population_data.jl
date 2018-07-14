#!/usr/bin/julia
#
# Parse county population data files
#
import JSON
include(dirname(@__FILE__) * "/common.jl")

# Get state codes
println("Parsing state codes...")
(state_codes, _) = readdlm(project_dir * "/data/state_codes.tsv", '\t', header=true)
#state_codes = [["California" "CA" 6]]
#state_codes = [["Texas" "TX" 48]]

# Import population and voting data
# Note: Each county GEOID maps to the tuple:
#   (population, Dem votes, GOP votes)
println("Parsing population and voting data...")
(vote_table, _) = readdlm(project_dir * "/data/county_pop_vote.tsv", '\t', header=true)
vote_data = Dict{Int64, Tuple{Int64, Int64, Int64}}()
for row in 1:size(vote_table)[1]
    geoid     = Int64(vote_table[row, 1])
    pop       = Int64(vote_table[row, 6])
    dem_votes = Int64(vote_table[row, 3])
    gop_votes = Int64(vote_table[row, 4])
    vote_data[geoid] = (pop, dem_votes, gop_votes)
end

# Parse census block data for each state
county_data = ["geoid" "pop" "dem_votes" "gop_votes" "pos_x" "pos_y" "pos_z" "pos_var"]
const deg_to_rad   = pi / 180
isdir(download_dir) || mkdir(download_dir) # Create directory if needed
for row in 1:size(state_codes)[1]
    state_name   = state_codes[row, 1]
    state_abbrev = state_codes[row, 2]
    state_id     = Int64(state_codes[row, 3])

    # Download shapefile and generate JSON file, if needed
    file_base = "tabblock2010_" * lpad(state_id, 2, "0") * "_pophu"
    json_file = download_dir * "/" * file_base * ".json"
    if !isfile(json_file)
        println("Downloading census block data for " * state_name * "...")
        census_url = "ftp://ftp2.census.gov/geo/tiger/TIGER2010BLKPOPHU/"
        zip_file = download_dir * "/" * file_base * ".zip"
        shape_file = download_dir * "/" * file_base * ".shp"
        download(census_url * file_base * ".zip", zip_file)
        run(`unzip $zip_file -d $download_dir`)
        run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
    end

    # Parse JSON file
    println("Parsing census block data for " * state_name * "...")
    sums = Dict{Int64, Vector{Float64}}()
    geoid_list = []
    for block_json in JSON.parsefile(json_file)["features"]

        # Initialize county data
        props_json = block_json["properties"]
        county_id = parse(Int64, props_json["COUNTYFP10"])
        pop = Int64(props_json["POP10"])
        geoid = state_id * 1000 + county_id
        if !haskey(sums, geoid)
            sums[geoid] = zeros(7)
            push!(geoid_list, geoid)
        end
        county_sums = sums[geoid]

        # Approximate census block position
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        num_coords = 0
        geometry = block_json["geometry"]["type"]
        if geometry == "Polygon"
            for polygon in block_json["geometry"]["coordinates"]
                for coord in polygon
                    long = deg_to_rad * Float64(coord[1])
                    lat  = deg_to_rad * Float64(coord[2])
                    sum_x += cos(long) * sin(lat)
                    sum_y += sin(long) * sin(lat)
                    sum_z += cos(lat)
                    num_coords += 1
                end
            end
        elseif geometry == "MultiPolygon"
            for multipolygon in block_json["geometry"]["coordinates"]
                for polygon in multipolygon
                    for coord in polygon
                        long = deg_to_rad * Float64(coord[1])
                        lat  = deg_to_rad * Float64(coord[2])
                        sum_x += cos(long) * sin(lat)
                        sum_y += sin(long) * sin(lat)
                        sum_z += cos(lat)
                        num_coords += 1
                    end
                end
            end            
        end
        x = earth_radius * sum_x / num_coords
        y = earth_radius * sum_y / num_coords
        z = earth_radius * sum_z / num_coords

        # Record population and position of census block
        county_sums[1] += pop * x
        county_sums[2] += pop * y
        county_sums[3] += pop * z
        county_sums[4] += pop * x * x
        county_sums[5] += pop * y * y
        county_sums[6] += pop * z * z
        county_sums[7] += pop
        
    end

    # Compute population mean and standard deviation
    state_county_data = Array{Any, 2}(length(geoid_list), 8)
    for (i, geoid) in enumerate(geoid_list)
        county_sums = sums[geoid]
        if !haskey(vote_data, geoid)
            println("GEOID " * string(geoid) * " not found in voting data")
        end
        county_votes = vote_data[geoid]
        pop = county_sums[7]
        x = county_sums[1] / pop
        y = county_sums[2] / pop
        z = county_sums[3] / pop
        var_x = county_sums[4] / pop - x * x
        var_y = county_sums[5] / pop - y * y
        var_z = county_sums[6] / pop - z * z
        var = max(var_x + var_y + var_z, 0.0)
        state_county_data[i, 1] = geoid
        state_county_data[i, 2] = county_votes[1]
        state_county_data[i, 3] = county_votes[2]
        state_county_data[i, 4] = county_votes[3]
        state_county_data[i, 5] = x
        state_county_data[i, 6] = y
        state_county_data[i, 7] = z
        state_county_data[i, 8] = var
    end
    county_data = [county_data; state_county_data]
    
end

# Output results to file
println("Writing results to file...")
isdir(results_dir) || mkdir(results_dir) # Create directory if needed
writedlm(county_data_file, county_data, '\t')
