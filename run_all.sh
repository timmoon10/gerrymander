#!/bin/bash
pushd $(dirname $(realpath $0))
#julia parse_population_data.jl
julia parse_geography_data.jl
julia construct_interaction_graph.jl
julia construct_geography_graph.jl
julia relax_partitions_softmax.jl $@
JULIA_NUM_THREADS=4 julia plot_results.jl
julia simulate_election.jl
popd
