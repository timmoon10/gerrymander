#!/usr/bin/julia
#
# Parse county geography data files
#
import JSON
using ProtoBuf
include(dirname(@__FILE__) * "/common.jl")
include(proto_file)

# Construct protobuf objects
println("Initializing protobuf objects...")
isdir(results_dir) || mkdir(results_dir)
if !isfile(results_dir * "/gerrymander_pb.jl")
    run(`$protoc --plugin=$julia_protobuf_dir/plugin/protoc-gen-julia 
         -I=$project_dir --julia_out=$results_dir
         gerrymander.proto`)
end
include(results_dir * "/gerrymander_pb.jl")
bounds_list_proto = CountyBoundariesList()
fillset(bounds_list_proto, :county_bounds)

# Download shapefile and generate JSON file, if needed
file_base = "tl_2010_us_county10"
json_file = download_dir * "/" * file_base * ".json"
if !isfile(json_file)
    println("Downloading county geography data...")
    census_url = "ftp://ftp2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/"
    zip_file = download_dir * "/" * file_base * ".zip"
    shape_file = download_dir * "/" * file_base * ".shp"
    download(census_url * file_base * ".zip", zip_file)
    run(`unzip $zip_file -d $download_dir`)
    run(`ogr2ogr -f "GeoJSON" $json_file $shape_file`)
end

# Parse GeoJSON data for each county
println("Parsing county geography data...")
for county_json in JSON.parsefile(json_file)["features"]
    bounds_proto = CountyBoundaries()
    geoid = parse(UInt32, county_json["properties"]["GEOID10"])
    set_field!(bounds_proto, :geoid, geoid)
    fillset(bounds_proto, :polygon)
    geometry_type = county_json["geometry"]["type"]
    if geometry_type == "Polygon"
        for polygon_coords in county_json["geometry"]["coordinates"]
            polygon_proto = Polygon()
            fillset(polygon_proto, :x)
            fillset(polygon_proto, :y)
            for coord in polygon_coords
                add_field!(polygon_proto, :x, coord[1])
                add_field!(polygon_proto, :y, coord[2])
            end
            add_field!(bounds_proto, :polygon, polygon_proto)
        end
    end
    if geometry_type == "MultiPolygon"
        for multipolygon in county_json["geometry"]["coordinates"]
            for polygon_coords in multipolygon
                polygon_proto = Polygon()
                fillset(polygon_proto, :x)
                fillset(polygon_proto, :y)
                for coord in polygon_coords
                    add_field!(polygon_proto, :x, coord[1])
                    add_field!(polygon_proto, :y, coord[2])
                end
                add_field!(bounds_proto, :polygon, polygon_proto)
            end
        end
    end
    add_field!(bounds_list_proto, :county_bounds, bounds_proto)
end

# Write results to file
println("Writing results to file...")
isdir(results_dir) || mkdir(results_dir)
open(county_bounds_file, "w") do f
    writeproto(f, bounds_list_proto)
end
