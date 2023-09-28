using Serialization
using Base.Threads
include("coevolution_network_base.jl")
using .CoevolutionNetworkBase

const DIRECTORY_PATH = "C:/Users/Daniel/Desktop/simresults_oneinfected/" # or wherever your files are located

function clean_data(directory_path::String)
    # Get a list of all files in the directory
    files = readdir(directory_path)
    
    Threads.@threads for file in files
        filepath = joinpath(directory_path, file) # Full path to the file
        
        # Try to open the file
        try
            open(deserialize, filepath)
        catch e
            println("Error processing file $filepath: ", e)
            
            # Delete the problematic file
            try
                rm(filepath)
                println("Deleted problematic file: $filepath")
            catch deletion_error
                println("Failed to delete file $filepath. Error: ", deletion_error)
            end
        end
    end
end

clean_data(DIRECTORY_PATH)
