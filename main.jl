using Random
using Printf
using Dates
include("./needleman_wunsch.jl")
Random.seed!(1)

debug = false

if debug
    function dprintf(s::String, args...)
        @eval @printf($s, $(args...))
    end
else
    function dprintf(s, args...)
        # nothing
    end
end

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
    dprintf("Comparing %s to %s\n", string(A), string(B))
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

function random_permutations(A::Array{String,1}, num::Int64)::Array{Array{String,1},1}
    # Return num random permutations of the elements of A
    results = Array{Array{String,1},1}(undef,0)
    if num >= factorial(big(length(A)))
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
function get_swap_sequence(A::Array{String,1}, B::Array{String,1})
    if length(A) != length(B)
        throw(ErrorException("Length of A and B must be the same"))
    end

    localA = Array{String,1}(undef,0)
    localB = Array{String,1}(undef,0)
    for i in 1:length(A)
        push!(localA, A[i])
        push!(localB, B[i])
    end

    # localA = deepcopy(A)
    # localB = deepcopy(B)
    swap_sequence = []

    for i in 1:length(A)
        dprintf("Checking index %d, localA: %s\nlocalB: %s\n", i, string(localA), string(localB))
        if localA[i] != localB[i]

            # check where A[i] is in B[i:end]
            b_indexes = indexin([localA[i]], localB[i:end])
            @assert(length(b_indexes) == 1)
            b_index = b_indexes[1]
            if b_index === nothing
                msg = string("Could not find ", string(localA[i], string(" in ", string(localB[i:end]))))
                throw(ErrorException(msg))
            end
            b_index += i - 1

            dprintf("Doing swap (%d,%d)\n", i, b_index)
            # record the swap of i and b_index
            temp = localA[i]
            localA[i] = localA[b_index]
            localA[b_index] = temp

            push!(swap_sequence, (i, b_index))
        end
    end
    dprintf("localA: %s\nlocalB: %s\n", string(localA), string(localB))
    local_equals = array_equals_ordered(localA, localB)
    dprintf("localA == localB: %s\n", string(local_equals))
    return swap_sequence
end

# Calculates a new "velocity" for the particle given global best particle, local best particle
function pso_particle_velocity(particle::Array{String,1}, global_best::Array{String,1}, local_best::Array{String,1}, alpha::Float64, beta::Float64)
    swaps_to_local_best = get_swap_sequence(particle, local_best)
    swaps_to_global_best = get_swap_sequence(particle, global_best)
    @printf("swaps_to_global_best: %s\n", swaps_to_global_best)
    @printf("swaps_to_local_best: %s\n", swaps_to_local_best)
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


function apply_velocity_to_particle(particle_in::Array{String,1}, swap_sequence_velocity::Array)::Array{String,1}
    particle = deepcopy(particle_in)
    for (idxA, idxB) in swap_sequence_velocity
        temp = particle[idxA]
        particle[idxA] = particle[idxB]
        particle[idxB] = temp
    end
    return particle
end

function get_new_gap_indexes(old::String, new::String)::Array{Int64,1}

    old_local = deepcopy(old)
    indexes = []
    for i in 1:length(new)
        if (new[i] == '-' && i == (length(old_local) + 1)) || new[i] == '-' && old_local[i] != '-'
            #@printf("get_new_gap_indexes, new gap at %d, old: %s, new: %s\n", i, old_local, new)
            push!(indexes, i)
            old_local = string(string(old_local[1:i-1], '-'), old_local[i:end])
        elseif new[i] == old_local[i]
            continue
        else
            #@printf("Sequences were not the same: %s, %s\n", old, new)
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
            #@printf("Inserting gap into %s at index %d\n", s, i)
            s = string(string(s[1:i-1], '-'), s[i:end])
            #@printf("New value: %s\n", s)
        end
        # put the modified string back in the array
        sequences[seq_index] = s
    end
end


function Make_Profile(k::Int64, t::Int64, A::Array{String,1})
    idxA = zeros(Int64, k)
    idxC = zeros(Int64, k)
    idxG = zeros(Int64, k)
    idxT = zeros(Int64, k)
    idxGap = zeros(Int64, k)

    for index in range(1, length=k)
        for seq in A
            if seq[index] == 'A'
                idxA[index] += 1
            elseif seq[index] == 'C'
                idxC[index] += 1
            elseif seq[index] == 'G'
                idxG[index] += 1
            elseif seq[index] == 'T'
                idxT[index] += 1
            elseif seq[index] == '-'
                idxGap[index] += 1
            else
                msg = @sprintf("character [%s] not allowed in DNA sequence", string(seq[index]))
                throw(ErrorException(msg))
            end
        end
    end

    return Dict(
        "A" => idxA,
        "C" => idxC,
        "G" => idxG,
        "T" => idxT,
        "-" => idxGap
    )
end

function score_sequences(A::Array{String,1})::Int64
    k = length(A[1])    # how long are the sequences
    t = length(A)       # how many sequences
    profile = Make_Profile(k, t, A)
    # score is the sum of the maximum occurring letter (not including gaps) in each position
    score = 0
    for i in 1:k
        score += max(profile["A"][i], profile["T"][i], profile["G"][i], profile["C"][i])
    end
    return score
end


function sequence_equals_ignore_gaps(A::String, B::String)
    A = filter(x -> x != '-', A)
    B = filter(x -> x != '-', B)
    return A == B
end

function progressive_alignment_inorder(sequences::Array, edges::Array{Tuple{String,String,Int},1})::Array{String,1}
    if length(edges) != length(sequences) - 1
        msg = @sprintf("Expected %d edge weights, got %d", length(sequences) - 1, length(edges))
        throw(ErrorException(msg))
    end
    local_edges = deepcopy(edges)
    total_score = 0

    aligned_sequences = Array{String,1}(undef,0)
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

        #@printf("aligned_sequences: %s\n", string(aligned_sequences))
        #@printf("original edges: %s\n", string(edges))
        #@printf("local_edges: %s\n", string(local_edges))
        dprintf("min_edge: %s\n", string(min_edge))

        # do the alignment of min_i and (min_i + 1)
        # TODO use aligned sequences here instead of original sequences
        #@printf("Looking for %s in %s\n", min_edge[1], sequences)
        idxA = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[1]), sequences)
        #idxA = findfirst(s -> s == min_edge[1], sequences)
        #@printf("Looking for %s in %s\n", min_edge[2], sequences)
        idxB = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[2]), sequences)
        #idxB = findfirst(s -> s == min_edge[2], sequences)
        A = aligned_sequences[idxA]
        B = aligned_sequences[idxB]
        @printf("\n\nPerforming global alignment of %d %d\n%s\nand\n%s\nedge_weight: %d\n", idxA, idxB, A, B, min_edge[3])
        score, alignedA, alignedB = global_align(A, B)
        @printf("Got alignment\n%s\nand\n%s\n", alignedA, alignedB)

        # TODO keep cumulative alignment.
        # need to figure out how to use alignedA and alignedB for the next step
        # TODO take new_gap_indexes_A, and insert same gaps into all current
        # predecessors of A.  And same for B
        new_gap_indexes_A = get_new_gap_indexes(A, alignedA)
        #@printf("new_gap_indexes_A: %s\n", string(new_gap_indexes_A))
        #@printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))
        seqs_to_update = copy(sets[idxA])
        deleteat!(seqs_to_update, findfirst(x -> x == idxA, seqs_to_update))
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_A)
        #@printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))

        new_gap_indexes_B = get_new_gap_indexes(B, alignedB)
        #@printf("new_gap_indexes_B: %s\n", string(new_gap_indexes_B))
        #println()
        #@printf("predecessors B(%d): %s\n", idxB, string(sets[idxB]))
        seqs_to_update = copy(sets[idxB])
        deleteat!(seqs_to_update, findfirst(x -> x == idxB, seqs_to_update))
        insert_gaps_at(aligned_sequences, sets[idxB], new_gap_indexes_B)
        #@printf("predecessors B(%d): %s\n", idxB, string(sets[idxB]))

        aligned_sequences[idxA] = alignedA
        aligned_sequences[idxB] = alignedB

        # for seq in aligned_sequences
        #     println(seq)
        # end

        # Update predecessors for A and B
        # TODO avoid the loop somehow by reusing the sets array from one of the existing idxA or idxB elements
        newset = []#[idxA, idxB]
        push!(newset, sets[idxA]...)
        push!(newset, sets[idxB]...)
        #@printf("newset: %s\n", string(newset))
        sets[idxA] = newset
        sets[idxB] = newset
        for p in newset
            sets[p] = newset
        end

    end
    return aligned_sequences
end

# Return the (edge, index)
function get_edge_between(nodeA::String, nodeB::String, edges::Array)
    nodeA = filter(x -> x != '-', nodeA)
    nodeB = filter(x -> x != '-', nodeB)

    for idx in 1:length(edges)
        edge = edges[idx]
        if edge[1] == nodeA && edge[2] == nodeB
            return (edge, idx)
        elseif edge[1] == nodeB && edge[2] == nodeA
            return (edge, idx)
        end
    end
    return nothing
end

# TODO maybe don't need SequenceMap
mutable struct SequenceMap
    original::String
    original_index::Int
    aligned::String
    aligned_index::Int
end

mutable struct Particle
    sequences::Array{String,1}
    aligned_sequences::Array{String,1}
    #sequence_map::Array{SequenceMap,1}
    #position::Array{Int,1} # ordering
    best_score::Int
    pidx::Int
    edges::Array{Tuple{String,String,Int},1}
end

function PSO_MSA()
    N = 40 # length of the generated sequences
    t = 10 # number of sequences
    iterations = 100
    solution_space = factorial(big(N)) / (N*(N-1))
    search_space = N * iterations
    search_solution_ratio = search_space / solution_space
    @printf("Solution Space: %.02f\n", solution_space)
    @printf("Search Space: %.02f\n", search_space)
    @printf("Search/Solution space: %.02f\n", search_solution_ratio)

    sequences = generate_sequences(t, N)
    println(sequences)

    nodes = Array{String,1}(undef,0)
    for i in 1:length(sequences)
        push!(nodes, sequences[i])
    end

    num_particles = 0
    if N*N > 100 || factorial(big(N*N)) > 100
        num_particles = 100
    else
        num_particles = N*N # num edges
    end

    permutations = random_permutations(nodes, t)

    # Initialize particles list from the generated permutations
    println("Initializing particles")
    particles = Array{Particle,1}(undef,length(permutations))
    for pidx in 1:length(permutations)
        sequences = copy(permutations[pidx])
        aligned_sequences = copy(sequences)
        # sequence_maps = Array{SequenceMap,1}(undef,length(sequences))
        # for sidx in 1:length(sequences)
        #     sequence_map = SequenceMap(sequences[sidx], sidx, sequences[sidx], sidx)
        #     sequence_maps[sidx] = sequence_map
        # end
        edges = Array{Tuple{String,String,Int},1}(undef,length(sequences)-1)
        initial_best = score_sequences(sequences)
        p = Particle(permutations[pidx], copy(permutations[pidx]), initial_best, pidx, edges)
        particles[pidx] = p
    end

    alpha = Random.rand()
    beta = Random.rand()

    # keep track of the local best score+ordering for each particle
    local_best_scores = [p.best_score for p in particles]
    local_best_particles = deepcopy(particles)
    # keep track of the global best particle
    global_best_particle = particles[argmax(local_best_scores)]
    @printf("Initial global best order:\n")
    for seq in global_best_particle.sequences
        println(seq)
    end
    @printf("Initial global best score: %d\n", global_best_particle.best_score)

    # Pick a starting position, arbitrary
    position_xid = copy(particles[1].sequences)

    iteration_metrics = []

    for iteration in 1:iterations
        # Revert all particle aligned sequences to their original unaligned values
        for p in particles
            p.aligned_sequences = copy(p.sequences)
        end

        @printf("Beginning iteration %d\n", iteration)
        start = Dates.now()

        # calculate the edges for each particle's sequences
        # only calculate the distance for the neighbors
        for pidx in 1:length(particles)
            particle = particles[pidx]
            p_edges = Array{Tuple{String,String,Int},1}(undef,0)
            for sidx in 2:length(particle.sequences)
                seqA = particle.sequences[sidx-1]
                seqB = particle.sequences[sidx]
                scoreAB = score_sequences([seqA, seqB])
                push!(p_edges, (seqA, seqB, scoreAB))
            end
            particle.edges = p_edges
        end

        # Filter list of particles to those only in same position as xid
        filtered_indexes = []
        for pidx in 1:length(particles)
            if array_equals_ordered(particles[pidx].sequences, position_xid)
                push!(filtered_indexes, pidx)
            end
        end

        # Get the new swap sequence for the particles at Xid.
        # Here position_xid is the current position as well as the local best,
        # since no particle can be in a position unless that position
        # is an improved score over the previous positions it has been in.
        new_velocity = pso_particle_velocity(
                position_xid,
                global_best_particle.sequences,
                position_xid, # this is the local best
                alpha, beta)

        if length(filtered_indexes) == 0
            @printf("filtered_indexes is empty, position_xid = %s\n", position_xid)
            continue
        else
            @printf("updating %d particles\n", length(filtered_indexes))
        end

        selected_particle_index = filtered_indexes[1]
        selected_particle = particles[selected_particle_index]

        dprintf("Performing progressive alignment of particle:\n%s\nEdges: %s\n", string(position_xid), string(selected_particle.edges))
        aligned_sequences = progressive_alignment_inorder(position_xid, selected_particle.edges)
        score = score_sequences(aligned_sequences)

        # apply the velocity to the TSP order

        new_position_xid = apply_velocity_to_particle(position_xid, new_velocity)
        @printf("Applied velocity %s to %s, and got %s\n", string(new_velocity), position_xid, new_position_xid)
        position_xid = new_position_xid

        # Update score for each particle in the current position, if the score is better
        for fidx in filtered_indexes
            p = particles[fidx]
            if score > p.best_score
                dprintf("New score for %s %d is better than previous score %d\n", position_xid, score, p.best_score)
                p.best_score = score
                p.aligned_sequences = copy(aligned_sequences)
                p.sequences = position_xid

                # local_best_scores[fidx] = score
                # local_best_particles[fidx] = position_xid
                # particles[fidx] = position_xid
            end
            if fidx == global_best_particle.pidx
                # update the aligned sequences
                global_best_particle.aligned_sequences = copy(aligned_sequences)
            end
        end

        if score > global_best_particle.best_score
            dprintf("New score for %s %d is better than previous global best %d\n",
                    position_xid, score, global_best_particle.best_score)
            global_best_particle = deepcopy(selected_particle)
            global_best_particle.best_score = score
            global_best_particle.aligned_sequences = copy(aligned_sequences)
        end


        end_time = Dates.now()
        duration = end_time - start
        push!(iteration_metrics, duration)
        @printf("\nIteration %d took %d milliseconds\n", iteration, duration.value)
        println("Global best:")
        for seq in global_best_particle.aligned_sequences
            println(seq)
        end
        println()
    end

    println("FINISHED PARTICLE SWARM")

    # print the final results
    @printf("\nGlobal best (score = %d):\n", global_best_particle.best_score)
    for seq in global_best_particle.aligned_sequences
        println(seq)
    end
    println()

    # remove any all-gap columns
    i = 1
    while i <= length(global_best_particle.aligned_sequences[1])
    #for i in 1:length(global_best_particle_value[1])
        #@printf("i = %d\n", i)
        all_gaps = true
        for seq in global_best_particle.aligned_sequences
            if seq[i] != '-'
                all_gaps = false
                break
            end
        end
        if all_gaps
            println("removing gap-only column")
            # remove the column at index i
            for seqidx in 1:length(global_best_particle.aligned_sequences)
                seq = global_best_particle.aligned_sequences[seqidx]
                seq = string(seq[1:i-1], seq[i+1:end])
                global_best_particle.aligned_sequences[seqidx] = seq
            end
            i -= 1
            #println(i)
        end
        i += 1
    end
    global_best_particle.best_score = score_sequences(global_best_particle.aligned_sequences)
    @printf("\nGlobal best AFTER REMOVING GAP COLUMNS (score = %d):\n", global_best_particle.best_score)
    for seq in global_best_particle.aligned_sequences
        println(seq)
    end
    println()

end

perf = @timed PSO_MSA()

println(perf)

# swap_sequence = get_swap_sequence([1,2,3,4,5,6], [4,3,2,1,6,5])
# println(swap_sequence)