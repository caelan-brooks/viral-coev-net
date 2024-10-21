all_adjacency_matrices = Vector{Matrix}(undef, 11)
all_adjacency_matrices[1] = [0 1 1 1; 
                             1 0 1 1;
                             1 1 0 1;
                             1 1 1 0]

all_adjacency_matrices[2] = [0 1 0 0; 
                             1 0 1 0;
                             0 1 0 1;
                             0 0 1 0]

all_adjacency_matrices[3] = [0 1 1 0; 
                             1 0 0 0;
                             1 0 0 1;
                             0 0 1 0]

all_adjacency_matrices[4] = [0 1 0 1; 
                             1 0 1 0;
                             0 1 0 1;
                             1 0 1 0]

all_adjacency_matrices[5] = [0 1 0 0; 
                             1 0 1 1;
                             0 1 0 1;
                             0 1 1 0]

all_adjacency_matrices[6] = [0 1 1 1; 
                             1 0 0 0;
                             1 0 0 1;
                             1 0 1 0]

all_adjacency_matrices[7] = [0 1 1 0; 
                             1 0 1 0;
                             1 1 0 1;
                             0 0 1 0]

all_adjacency_matrices[8] = [0 1 0 0; 
                             1 0 1 1;
                             0 1 0 0;
                             0 1 0 0]

all_adjacency_matrices[9] = [0 1 1 1; 
                             1 0 0 0;
                             1 0 0 0;
                             1 0 0 0]

all_adjacency_matrices[10]= [0 1 1 1; 
                             1 0 1 0;
                             1 1 0 1;
                             1 0 1 0]

all_adjacency_matrices[11]= [0 1 1 0; 
                             1 0 1 1;
                             1 1 0 1;
                             0 1 1 0]

# Function to check if a graph is connected
function is_connected(matrix)
    M2 = matrix * matrix
    M3 = M2 * matrix
    sum_matrix = matrix + M2 + M3
    return all(sum_matrix .> 0)  # Checks if all entries are greater than 0
end

# Function to check if two matrices are the same
function matrices_are_same(mat1, mat2)
    return all(mat1 .== mat2)
end

# Check that no two matrices are the same
for i in 1:length(all_adjacency_matrices)
    for j in i+1:length(all_adjacency_matrices)
        if matrices_are_same(all_adjacency_matrices[i], all_adjacency_matrices[j])
            println("Matrix $i and Matrix $j are the same.")
        end
    end
end

# Check if each graph is connected
for (i, matrix) in enumerate(all_adjacency_matrices)
    if is_connected(matrix)
        println("Graph $i is connected.")
    else
        println("Graph $i is not connected.")
    end
end

function normalize_matrix(matrix)
    column_sums = sum(matrix, dims=1)  # Sums across rows to get the column sums
    return matrix ./ column_sums
end

# Function to check if a matrix is symmetric
function is_symmetric(matrix)
    return matrix == matrix'
end

# Check symmetry, connectivity, and then normalize matrices
for i in 1:length(all_adjacency_matrices)
    adj_matrix = all_adjacency_matrices[i]
    if !is_symmetric(adj_matrix)
        println("Matrix $i is not symmetric.")
        continue  # Skip to the next matrix if this one is not symmetric
    end
    if !is_connected(adj_matrix)
        println("Graph $i (Matrix $i) is not connected.")
        continue  # Skip to the next matrix if this one is not connected
    end
    
    # Normalize the matrix
    # all_adjacency_matrices[i] = normalize_matrix(adj_matrix)
    # println("Matrix $i normalized.")

    all_adjacency_matrices[i] = adj_matrix
    println("Matrix $i is not normalized")
end