# Parse county election returns in 2020 presidential election
# Note: This requires downloading "County Presidential Election
# Returns 2000-2020" from MIT Election Data and Science Lab
# (https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ).

import DataStructures
import DelimitedFiles

function main()

    # Load election results file
    if length(ARGS) < 1
        throw("Must provide path to election results file")
    end
    in_file = ARGS[1]
    if !isfile(in_file)
        throw("Could not file election results file at $data_file")
    end
    (raw_data, _) = DelimitedFiles.readdlm(in_file, ',', header=true)

    # Parse election results file
    make_dict = () -> Dict{String, Int}("DEMOCRAT" => 0, "REPUBLICAN" => 0, "OTHER" => 0)
    votes = DataStructures.DefaultDict{Int, Dict{String, Int}}(make_dict)
    for i in 1:size(raw_data, 1)
        if Int(raw_data[i,1]) != 2020
            continue
        end
        state = String(raw_data[i,3])
        if state == "AK" || state == "HI"
            # Exclude non-contiguous states
            continue
        end
        county::Int = 0
        try
            county = Int(raw_data[i,5])
        catch err
            # Some precincts do not have FIPS code
            continue
        end
        party = String(raw_data[i,8])
        if party != "DEMOCRAT" && party != "REPUBLICAN"
            party = "OTHER"
        end
        votes[county][party] += Int(raw_data[i,9])
    end

    # Convert election results to table
    data = ["geoid" "dem" "gop" "other"]
    counties = collect(keys(votes))
    sort!(counties)
    for county in counties
        county_votes = votes[county]
        dem = county_votes["DEMOCRAT"]
        gop = county_votes["REPUBLICAN"]
        other = county_votes["OTHER"]
        data = [data; county dem gop other]
    end

    # Output results to file
    root_path = realpath(dirname(realpath(@__DIR__)))
    out_file = joinpath(root_path, "data", "election2020.csv")
    DelimitedFiles.writedlm(out_file, data, ',')

end

# Main function
main()
