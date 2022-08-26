#!/bin/bash
pushd $(dirname $(realpath $0))
julia relax_partitions_softmax.jl $@
JULIA_NUM_THREADS=4 julia plot_results.jl
popd
