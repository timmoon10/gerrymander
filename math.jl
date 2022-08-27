#
# Math utility functions
#

function softmax!(x::AbstractArray{Float64,1})
    shift = -maximum(x)
    for i in 1:length(x)
        x[i] = exp(x[i] + shift)
    end
    scale = 1 / sum(x)
    for i in 1:length(x)
        x[i] *= scale
    end
end

function zscore!(x::AbstractArray{Float64,1})
    mean::Float64 = sum(x) / length(x)
    var::Float64 = 0
    for i in 1:length(x)
        diff::Float64 = x[i] - mean
        var += diff * diff
    end
    var /= length(x)
    rstd::Float64 = 1 / sqrt(var)
    x .-= mean
    x .*= rstd
end
