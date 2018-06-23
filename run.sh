#!/bin/sh
project_dir=/home/moon/Documents/gerrymander/
rm -r $project_dir/output
$project_dir/parse_data.jl
$project_dir/construct_graph.jl
$project_dir/partition_counties.jl
$project_dir/plot_results.jl
