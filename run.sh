#!/bin/bash
pushd $(dirname $(realpath $0))
julia relax_partitions_softmax.jl $@
julia plot_results.jl
popd
