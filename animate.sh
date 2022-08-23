#!/bin/bash

# Options
num_steps=10
interval=1
num_partitions=6

# Apply partition relaxation steps
pushd $(dirname $(realpath $0))
rm results/partition.tsv
rm -r results/animation
mkdir -p results/animation
for (( step = 1; step <= ${num_steps}; step += ${interval} )); do
    julia \
        relax_partitions.jl \
        --relaxation-steps=${interval} \
        --num-partitions=${num_partitions}
    julia plot_results.jl
    cp results/partition.png results/animation/${step}.png
done
popd
