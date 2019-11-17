using Random
using Printf
include("./needleman_wunsch.jl")
Random.seed!(1)

function generate_sequences(t::Int64, l::Int64)
    # t is the number of sequences to create
    # l is the length of the sequences
    DNA = Array{String,1}(undef,0)
    base_arr = ["A", "T", "G", "C"]

    for t_index in 1:t
        push!(DNA, "")
        for l_value in 1:l
            r = convert(Int64, floor(Random.rand() * 4) + 1)
            DNA[t_index] = string(DNA[t_index], base_arr[r])
        end
    end
    return DNA
end


function array_equals_ordered(A::Array, B::Array)
    @printf("Comparing %s to %s\n", string(A), string(B))
    if length(A) != length(B)
        return false
    end

    for i in 1:length(A)
        if A[i] != B[i]
            return false
        end
    end
    return true
end

function array_contains(A::Array, e)
    for i in 1:length(A)
        if array_equals_ordered(A[i], e)
            return true
        end
    end

    return false
end

function random_permutations(A, num::Int64)
    # Return num random permutations of the elements of the Array A
    results = []
    if num >= factorial(length(A))
        throw(ErrorException("num is too high"))
    end
    max_check = 50
    try_count = 0
    for i in 1:num
        #@printf("Generating permutation %d\n", i)
        #@printf("Current Results: %s\n", string(results))
        localA = deepcopy(A)
        randomA = shuffle!(localA)
        #@printf("Sequence: %s\n", string(randomA))
        while array_contains(results, randomA)
            randomA = shuffle!(A)
            #@printf("Sequence: %s\n", string(randomA))
            try_count += 1
            if try_count > max_check
                throw(ErrorException("Exceeded max permutation tries"))
            end
        end
        push!(results, randomA)
    end
    return results
end

# Returns a swap sequence [(idx1, idx2), ...] to turn A into B
# Not necessarily the optimal swap sequence
function get_swap_sequence(A::Array, B::Array)
    if length(A) != length(B)
        throw(ErrorException("Length of A and B must be the same"))
    end

    localA = deepcopy(A)
    localB = deepcopy(B)
    swap_sequence = []

    for i in 1:length(A)
        @printf("Checking index %d\n", i)
        @printf("localA: %s\nlocalB: %s\n", string(localA), string(localB))
        if localA[i] != localB[i]

            # check where A[i] is in B[i:end]
            b_indexes = indexin(localA[i], localB[i:end])
            @assert(length(b_indexes) == 1)
            b_index = b_indexes[1]
            if b_index === nothing
                msg = string("Could not find ", string(localA[i], string(" in ", string(localB[i:end]))))
                throw(ErrorException(msg))
            end
            b_index += i - 1

            @printf("Doing swap (%d,%d)\n", i, b_index)
            # record the swap of i and b_index
            temp = localA[i]
            localA[i] = localA[b_index]
            localA[b_index] = temp

            push!(swap_sequence, (i, b_index))
        end
    end
    @printf("localA: %s\nlocalB: %s\n", string(localA), string(localB))
    local_equals = array_equals_ordered(localA, localB)
    @printf("localA == localB: %s\n", string(local_equals))
    return swap_sequence
end

# Calculates a new "velocity" for the particle given global best score, local best score
function pso_particle_velocity(particle, global_best::Array, local_best::Array, alpha::Float64, beta::Float64)
    swaps_to_local_best = get_swap_sequence(particle, local_best)
    swaps_to_global_best = get_swap_sequence(particle, global_best)
    r = Random.rand()
    velocity = []
    if r < alpha
        append!(velocity, swaps_to_local_best)
    end
    if r < beta
        append!(velocity, swaps_to_global_best)
    end
    return velocity
end

function get_new_gap_indexes(old::String, new::String)::Array{Int64,1}

    old_local = deepcopy(old)
    indexes = []
    for i in 1:length(new)
        if (new[i] == '-' && i == (length(old_local) + 1)) || new[i] == '-' && old_local[i] != '-'
            @printf("get_new_gap_indexes, new gap at %d, old: %s, new: %s\n", i, old_local, new)
            push!(indexes, i)
            old_local = string(string(old_local[1:i-1], '-'), old_local[i:end])
        elseif new[i] == old_local[i]
            continue
        else
            @printf("Sequences were not the same: %s, %s\n", old, new)
            throw(ErrorException("Sequences ignoring gaps were not the same"))
        end
    end
    return indexes
end

function insert_gaps_at(sequences::Array, seq_indexes::Array{Int64,1}, gap_indexes::Array{Int64,1})
    gap_indexes = sort(gap_indexes)

    for sidx in 1:length(seq_indexes)
        seq_index = seq_indexes[sidx]
        s = sequences[seq_index]

        for i in gap_indexes
            # insert a '-' character at index i
            @printf("Inserting gap into %s at index %d\n", s, i)
            s = string(string(s[1:i-1], '-'), s[i:end])
            @printf("New value: %s\n", s)
        end
        # put the modified string back in the array
        sequences[seq_index] = s
    end
end

function progressive_alignment_inorder(sequences::Array, edges::Array{Tuple,1})
    if length(edges) != length(sequences) - 1
        msg = @sprintf("Expected %d edge weights, got %d", length(sequences) - 1, length(edges))
        throw(ErrorException(msg))
    end
    local_edges = deepcopy(edges)
    total_score = 0

    aligned_sequences = []
    predecessors = Array{Array{Int64,1}}(undef,0) # store which sequence indexes are at or below each other sequence in the tree

    sets = Array{Array{Int64,1},1}(undef,0)

    for i in 1:length(sequences)
        #seq in sequences
        seq = sequences[i]
        push!(aligned_sequences, deepcopy(seq))
        push!(predecessors, [])

        push!(sets, [i])
    end


    for i in 1:length(edges)
        # find the min edge weight that has not been used yet
        min_edge = local_edges[1]
        min_edge_i = 1
        for edge_i in 2:length(local_edges)
            if local_edges[edge_i][3] < min_edge[3]
                min_edge = local_edges[edge_i]
                min_edge_i = edge_i
            end
        end
        deleteat!(local_edges, min_edge_i)

        # min_edge_weight = min(local_edge_weights...)
        # min_edge_i = Nothing
        # for orig_i in 1:length(edge_weights)
        #     if edge_weights[orig_i] == min_edge_weight
        #         min_edge_i = orig_i
        #         break
        #     end
        # end
        # edge_weight = edge_weights[min_edge_i]


        # #delete!(local_edge_weights, min_edge_weight)
        # deleteat!(local_edge_weights, findfirst(x -> x == min_edge_weight, local_edge_weights))

        #min_i += (i - 1)

        @printf("aligned_sequences: %s\n", string(aligned_sequences))
        @printf("original edges: %s\n", string(edges))
        @printf("local_edges: %s\n", string(local_edges))

        # do the alignment of min_i and (min_i + 1)
        idxA = findfirst(s -> s == min_edge[1], sequences)
        idxB = findfirst(s -> s == min_edge[2], sequences)
        A = aligned_sequences[idxA]
        B = aligned_sequences[idxB]
        @printf("\n\nPerforming global alignment of %s and %s, edge_weight: %d\n", A, B, min_edge[3])
        score, alignedA, alignedB = global_align(A, B)


        # TODO keep cumulative alignment.
        # need to figure out how to use alignedA and alignedB for the next step
        # TODO take new_gap_indexes_A, and insert same gaps into all current
        # predecessors of A.  And same for B
        new_gap_indexes_A = get_new_gap_indexes(A, alignedA)
        @printf("new_gap_indexes_A: %s\n", string(new_gap_indexes_A))
        @printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))
        seqs_to_update = copy(sets[idxA])
        deleteat!(seqs_to_update, findfirst(x -> x == idxA, seqs_to_update))
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_A)
        @printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))

        new_gap_indexes_B = get_new_gap_indexes(B, alignedB)
        @printf("new_gap_indexes_B: %s\n", string(new_gap_indexes_B))
        @printf("predecessors B(%d): %s\n", idxB, string(sets[idxB]))
        seqs_to_update = copy(sets[idxB])
        deleteat!(seqs_to_update, findfirst(x -> x == idxB, seqs_to_update))
        insert_gaps_at(aligned_sequences, sets[idxB], new_gap_indexes_B)
        @printf("predecessors B(%d): %s\n", idxB, string(sets[idxB]))

        aligned_sequences[idxA] = alignedA
        aligned_sequences[idxB] = alignedB

        # Update predecessors for A and B
        #curr_predecessors_A = deepcopy(predecessors[idxA])

        # TODO avoid the loop somehow by reusing the sets array from one of the existing idxA or idxB elements
        newset = []#[idxA, idxB]
        push!(newset, sets[idxA]...)
        push!(newset, sets[idxB]...)
        @printf("newset: %s\n", string(newset))
        sets[idxA] = newset
        sets[idxB] = newset
        for p in newset
            sets[p] = newset
        end

    end
    return aligned_sequences
end

function get_edge_between(nodeA::String, nodeB::String, edges::Array)
    for edge in edges
        if edge[1] == nodeA && edge[2] == nodeB
            return edge
        elseif edge[1] == nodeB && edge[2] == nodeA
            return edge
        end
    end
    return Nothing
end

function PSO_MSA()
    N = 10
    t = 5
    iterations = 10
    solution_space = factorial(N) / (N*(N-1))
    search_space = N * iterations
    search_solution_ratio = search_space / solution_space
    @printf("Solution Space: %.02f\n", solution_space)
    @printf("Search Space: %.02f\n", search_space)
    @printf("Search/Solution space: %.02f\n", search_solution_ratio)

    sequences = generate_sequences(t, N)
    println(sequences)

    nodes = []
    for i in 1:length(sequences)
        push!(nodes, sequences[i])
    end

    edges = []
    for i in 1:length(sequences)
        for j in 1:length(sequences)
            if i != j
                seqA = sequences[i]
                seqB = sequences[j]
                # Use global_align function from needleman_wunsch to compute the distance
                distanceAB, seqA_aligned, seqB_aligned = global_align(seqA, seqB, 1, -1, 0)
                # Add edge A-B = distanceAB
                push!(edges, (seqA, seqB, distanceAB))
            end
        end
    end

    print("Nodes: ")
    println(nodes)
    print("Edges: ")
    println(edges)
    num_particles = 0
    if length(edges) > 100 || factorial(length(edges)) > 100
        num_particles = 100
    else
        num_particles = length(edges)
    end

    particles = random_permutations(nodes, 10)
    println("Permutations:")
    println(particles)

    alpha = Random.rand()
    beta = Random.rand()
    # TODO
    #global_best =  # calculate best permutation from the current particles
    #local_bests = zeros(Float64, length(particles))
    particle_scores = zeros(Float64, length(particles))

    #for iteration in 1:iterations

    # TODO get the edges between this particle ordering
    # look at particles[0] to test
    particles_1_edges = Array{Tuple,1}(undef,0)
    for i in 2:length(particles[1])
        edge = get_edge_between(particles[1][i-1], particles[1][i], edges)
        push!(particles_1_edges, edge) # push the edge weight only
    end
    @printf("Performing progressive alignment of particle 1:\n%s\nEdges: %s\n", string(particles[1]), string(particles_1_edges))
    aligned_sequences = progressive_alignment_inorder(particles[1], particles_1_edges)

    # TODO need function to find the "SCORE" of the aligned_sequences array

    # in each iteration, find the maximum score (minimize distance)

    for aligned_seq in aligned_sequences
        @printf("%s (length=%d)\n", aligned_seq, length(aligned_seq))
    end

end

PSO_MSA()

# swap_sequence = get_swap_sequence([1,2,3,4,5,6], [4,3,2,1,6,5])
# println(swap_sequence)