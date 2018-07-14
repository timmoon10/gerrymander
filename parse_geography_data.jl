#!/usr/bin/julia
#
# Parse geography data files
#
import JSON
using ProtoBuf
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Download shapefile and generate JSON file, if needed
file_base = "tl_2010_us_county10"
json_file = download_dir * "/" * file_base * ".json"
if !isfile(json_file)
    println("Downloading geography data...")
    census_url = "ftp://ftp2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/"
    zip_file = download_dir * "/" * file_base * ".zip"
    shape_file = download_dir * "/" * file_base * ".shp"
    download(census_url * file_base * ".zip", zip_file)
    run(`unzip $zip_file -d $download_dir`)
    run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
end

# Function to convert coordinates to a protobuf polygon
function coords_to_polygon_proto(polygon_coords)
    polygon_proto = Polygon()
    if length(polygon_coords) > 1
        fillset(polygon_proto, :interior_border)
    end
    for (i, border_coords) in enumerate(polygon_coords)
        border_proto = Ring()
        fillset(border_proto, :x)
        fillset(border_proto, :y)
        for coord in border_coords
            add_field!(border_proto, :x, coord[1])
            add_field!(border_proto, :y, coord[2])
        end
        if i == 1
            set_field!(polygon_proto, :exterior_border, border_proto)
        else
            add_field!(polygon_proto, :interior_border, border_proto)
        end
    end
    return polygon_proto
end

# Parse GeoJSON data
println("Parsing geography data...")
region_list_proto = MultiPolygonList()
fillset(region_list_proto, :multi_polygon)
for region_json in JSON.parsefile(json_file)["features"]
    id = parse(Int64, region_json["properties"]["GEOID10"])
    geometry_type = region_json["geometry"]["type"]
    region_proto = MultiPolygon()
    set_field!(region_proto, :id, id)
    fillset(region_proto, :polygon)
    if geometry_type == "Polygon"
        polygon_coords = region_json["geometry"]["coordinates"]
        polygon_proto = coords_to_polygon_proto(polygon_coords)
        add_field!(region_proto, :polygon, polygon_proto)
    elseif geometry_type == "MultiPolygon"
        for polygon_coords in region_json["geometry"]["coordinates"]
            polygon_proto = coords_to_polygon_proto(polygon_coords)
            add_field!(region_proto, :polygon, polygon_proto)
        end
    end
    add_field!(region_list_proto, :multi_polygon, region_proto)
end

# Write results to file
println("Exporting geography data...")
open(geography_data_file, "w") do f
    writeproto(f, region_list_proto)
end
