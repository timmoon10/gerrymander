#!/bin/sh
project_dir=$(git rev-parse --show-toplevel)
$project_dir/parse_population_data.jl
$project_dir/construct_interaction_graph.jl
$project_dir/parse_geography_data.jl
$project_dir/construct_geography_graph.jl
$project_dir/partition_counties.jl
$project_dir/plot_results.jl
