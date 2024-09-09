using Base.Threads
println("Number of threads: ", nthreads())
# Function to delete all files in a given directory
function delete_files_in_directory(dir)
    println(length(dir))
    flush(stdout)
    @threads for file in readdir(dir)
        full_path = joinpath(dir, file)
        if isfile(full_path)
            rm(full_path)
        elseif isdir(full_path)
            delete_files_in_directory(full_path)
        end
    end
    rm(dir)  # remove the subdirectory after clearing its contents
end

# Main function to handle multithreaded deletion
function delete_all_files_multithreaded(base_dir)
    # Collect all subdirectories two levels deep
    dirs_to_clean = []
    for subdir in readdir(base_dir)
        full_path = joinpath(base_dir, subdir)
        if isdir(full_path)
            for subsubdir in readdir(full_path)
                push!(dirs_to_clean, joinpath(full_path, subsubdir))
            end
        end
    end

    # Perform multithreaded deletion
    for dir in dirs_to_clean
        delete_files_in_directory(dir)
    end
end

# Base directory to clean
base_dir = "/pool001/dswartz/peak_scaling"

# Execute the cleanup
delete_all_files_multithreaded(base_dir)
