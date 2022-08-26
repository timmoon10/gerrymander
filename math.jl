#
# Math utility functions
#

function softmax!(x::Array{Float64})
    shift = -maximum(x)
    for i in 1:length(x)
        x[i] = exp(x[i] + shift)
    end
    scale = 1 / sum(x)
    for i in 1:length(x)
        x[i] *= scale
    end
end

function zscore!(x::Array{Float64})
    mean::Float64 = Statistics.mean(x)
    std::Float64 = Statistics.stdm(x, mean)
    rstd::Float64 = 1 / std
    x .-= mean
    x .*= rstd
end
