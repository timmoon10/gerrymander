module DataFiles

import DataStructures
import DelimitedFiles
import Memoize

export root_dir
@Memoize.memoize function root_dir()::String
    return realpath(dirname(realpath(@__DIR__)))
end

export state_names
@Memoize.memoize function state_names()::DataStructures.OrderedDict{UInt, String}
    (data, _) = DelimitedFiles.readdlm(
        joinpath(root_dir(), "data", "state_codes.tsv"),
        '\t',
        header=true,
        comments=true,
    )
    names = DataStructures.OrderedDict{UInt, String}()
    for row in 1:size(data, 1)
        names[data[row, 3]] = data[row, 1]
    end
    return names
end

end  # module DataFiles
