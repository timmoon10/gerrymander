#!/usr/bin/julia

#
# Parse data files for county geography, population, and votes
#
# TIGER/Line Shapefiles were obtained from the U.S. Census Bureau
# website (https://www.census.gov/geo/maps-data/data/tiger-line.html)
# and converted to GeoJSON format using GDAL (see
# https://stackoverflow.com/questions/2223979/convert-a-shapefile-shp-to-xml-json).
#
# Population data was obtained from
# https://en.wikipedia.org/wiki/List_of_United_States_counties_and_county_equivalents
# and 2016 voting data from
# https://github.com/tonmcg/County_Level_Election_Results_12-16.
#

# Parameters
desc = "California"
project_dir   = "/home/moon/Documents/gerrymander/"
output_dir    = project_dir * "/output"
shape_file    = "/home/moon/Downloads/tl_2017_us_county/us_county/tigerline.json"
pop_data_file = "/home/moon/Downloads/tl_2017_us_county/county_pops_votes.csv"
props_file    = output_dir * "/props.prototxt"  # Output file with county properties
bounds_file   = output_dir * "/bounds.prototxt" # Output file with county boundaries

# Import packages
import JSON
using ProtoBuf

# Import population and voting data
# Note: Each county GEOID maps to the tuple:
#   (population, Dem votes, GOP votes)
println("Parsing population and voting data...")
pop_data = Dict{UInt32, Tuple{UInt32, UInt32, UInt32}}()
pop_csv = readcsv(pop_data_file, skipstart=1)
for row in 1:size(pop_csv)[1]
    geoid = UInt32(pop_csv[row, 1])
    pop = UInt32(pop_csv[row, 6])
    dem_votes = UInt32(pop_csv[row, 3])
    gop_votes = UInt32(pop_csv[row, 4])
    pop_data[geoid] = (pop, dem_votes, gop_votes)
end

# Construct protobuf objects
println("Initializing protobuf objects...")
if !isfile(output_dir * "/gerrymander_pb.jl")
    protoc = "/home/moon/src/protobuf/src/protoc"
    julia_protobuf = "/home/moon/.julia/v0.4/ProtoBuf"
    run(`mkdir -p $output_dir`)
    run(`$protoc --plugin=$julia_protobuf/plugin/protoc-gen-julia 
         -I=$project_dir --julia_out=$output_dir
         gerrymander.proto`)
end
include(output_dir * "/gerrymander_pb.jl")
props_list_proto = CountyPropertiesList()
bounds_list_proto = CountyBoundariesList()
set_field!(props_list_proto, :desc, desc)
set_field!(bounds_list_proto, :desc, desc)
fillset(props_list_proto, :county_props)
fillset(bounds_list_proto, :county_bounds)

# Parse GeoJSON data for each county
println("Parsing geography data...")
for county_json in JSON.parsefile(shape_file)["features"]
    props_json = county_json["properties"]
    name = props_json["NAME"]
    geoid = parse(UInt32, props_json["GEOID"])
    stateid = parse(UInt32, props_json["STATEFP"])

    # Skip if county is not in continental U.S.
    if (stateid == 2  || stateid == 3  || stateid == 7  || stateid == 14 ||
        stateid == 15 || stateid == 43 || stateid == 52 || stateid >  56 )
        continue
    end
    
    # Convert county properties to protobuf
    props_proto = CountyProperties()
    county_pop_data = pop_data[geoid]
    set_field!(props_proto, :name, name)
    set_field!(props_proto, :geoid, geoid)
    set_field!(props_proto, :state_fips, stateid)
    set_field!(props_proto, :county_fips, parse(UInt32, props_json["COUNTYFP"]))
    set_field!(props_proto, :area, Float32(props_json["ALAND"] + props_json["AWATER"]))
    interior_point = Coord()
    set_field!(interior_point, :long, parse(Float32, props_json["INTPTLON"]))
    set_field!(interior_point, :lat, parse(Float32, props_json["INTPTLAT"]))
    set_field!(props_proto, :interior_point, interior_point)
    set_field!(props_proto, :population, county_pop_data[1])
    set_field!(props_proto, :dem_votes, county_pop_data[2])
    set_field!(props_proto, :gop_votes, county_pop_data[3])
    add_field!(props_list_proto, :county_props, props_proto)

    # Convert county boundaries to protobuf
    bounds_proto = CountyBoundaries()
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

# Write results to protobuf
println("Writing results to file...")
open(props_file, "w") do f
    writeproto(f, props_list_proto)
end
open(bounds_file, "w") do f
    writeproto(f, bounds_list_proto)
end
