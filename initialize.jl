#!/usr/bin/julia
#
# Initialize project directory
#
include(dirname(@__FILE__) * "/common.jl")

# Move to project directory
original_dir = pwd()
cd(project_dir)

# Create directories if needed
println("Creating directories...")
if !isdir(download_dir)
    mkdir(download_dir)
end
if !isdir(results_dir)
    mkdir(results_dir)
end

# Compile proto file
println("Initializing protobuf...")
julia_protobuf_plugin = Pkg.dir("ProtoBuf") * "/plugin/protoc-gen-julia"
proto_file_dir = dirname(proto_file)
run(`$protoc --plugin=$julia_protobuf_plugin -I=$project_dir
     --julia_out=$proto_file_dir $project_dir/gerrymander.proto`)

# Download submodules
println("Initializing Git submodules...")
run(`git submodule update --init --recursive`)

# Return to original directory
cd(original_dir)
