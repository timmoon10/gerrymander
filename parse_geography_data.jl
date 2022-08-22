#
# Parse geography data files
#
import DelimitedFiles
import Downloads
import JSON
import Serialization
include(joinpath(dirname(@__FILE__), "common.jl"))

# Import county data
println("Importing county data...")
(county_data, _) = DelimitedFiles.readdlm(
    county_data_file,
    '\t',
    header=true,
)
id_set = Set{Int64}()
for row in 1:size(county_data, 1)
    push!(id_set, county_data[row, 1])
end

# Download shapefile and generate JSON file, if needed
file_base = "tl_2010_us_county10"
json_file = joinpath(download_dir, file_base * ".json")
if !isfile(json_file)
    println("Downloading geography data...")
    census_url = "https://www2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/"
    zip_file = joinpath(download_dir, file_base * ".zip")
    shape_file = joinpath(download_dir, file_base * ".shp")
    Downloads.download(census_url * file_base * ".zip", zip_file)
    run(`unzip $zip_file -d $download_dir`)
    run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
end

# Data structure for geography data
# Note: Each dict key corresponds to a county. Each county is composed
# of one or more polygons, which have one or more borders. Each border
# is a list of (x,y) coordinates that form a complete ring. The first
# border is the external border, and other borders are internal.
regions = Dict{Int64, Vector{Vector{Array{Float64, 2}}}}()

# Function to parse JSON polygon
function parse_polygon(polygon_json::Vector{Any})::Vector{Array{Float64, 2}}
    polygon = Vector{Array{Float64, 2}}()
    for border_json::Vector{Any} in polygon_json
        border = Array{Float64, 2}(undef, (2, length(border_json)))
        for (i, coord::Vector{Float64}) in enumerate(border_json)
            border[1, i] = coord[1]
            border[2, i] = coord[2]
        end
        push!(polygon, border)
    end
    return polygon
end

# Parse GeoJSON data
println("Parsing geography data...")
for region_json in JSON.parsefile(json_file)["features"]
    id = parse(Int64, region_json["properties"]["GEOID10"])
    if !in(id, id_set)
        continue
    end
    geometry_type = region_json["geometry"]["type"]
    regions[id] = Vector{Vector{Array{Float64, 2}}}()
    if geometry_type == "Polygon"
        polygon_json = region_json["geometry"]["coordinates"]
        regions[id] = [parse_polygon(polygon_json)]
    elseif geometry_type == "MultiPolygon"
        regions[id] = []
        for polygon_json in region_json["geometry"]["coordinates"]
            push!(regions[id], parse_polygon(polygon_json))
        end
    end
end

# Write results to file
println("Exporting geography data...")
Serialization.serialize(geography_data_file, regions)
