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
    #dprintf("Comparing %s to %s\n", string(A), string(B))
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
        localA = deepcopy(A)
        randomA = shuffle!(localA)
        while array_contains(results, randomA)
            randomA = shuffle!(A)
            try_count += 1
            if try_count > max_check
                msg = @sprintf("Exceeded max permutation tries: %d", max_check)
                throw(ErrorException(msg))
            end
        end
        push!(results, randomA)
    end
    return results
end

# Returns a swap sequence [(idx1, idx2), ...] to turn A into B
# Not necessarily the optimal swap sequence
function get_swap_sequence(A_in::Array{String,1}, B_in::Array{String,1})
    if length(A_in) != length(B_in)
        throw(ErrorException("Length of A and B must be the same"))
    end
    #@printf("Getting swaps between:\n%s,\n%s\n", string(A_in), string(B_in))
    A = Array{String,1}(undef,0)
    B = Array{String,1}(undef,0)
    for i in 1:length(A_in)
        push!(A, A_in[i])
        push!(B, B_in[i])
    end

    swap_sequence = []

    for i in 1:length(A)
        #@printf("\nChecking index %d,\nlocalA: %s\nlocalB: %s\n", i, string(A), string(B))
        if A[i] != B[i]
            # We want to turn A into B
            # so find where B[i] is in A[i:end], and swap that element in A with the element at A[i]
            a_indexes = indexin([B[i]], A[i:end])
            @assert(length(a_indexes) == 1)
            a_index = a_indexes[1]
            if a_index === nothing
                msg = string("Could not find ", string(B[i], string(" in ", string(A[i:end]))))
                throw(ErrorException(msg))
            end

            a_index += i - 1

            #dprintf("Doing swap in A (%d, %d)\n", i, a_index)
            temp = A[i]
            A[i] = A[a_index]
            A[a_index] = temp

            push!(swap_sequence, (i, a_index))
        end
    end
    equals_check = array_equals_ordered(A, B)
    if !equals_check
        throw(ErrorException("A and B were not the same:\n%s\n%s\n", string(A), string(B)))
    end
    return swap_sequence
end

# Calculates a new "velocity" for the particle given global best particle, local best particle
function pso_particle_velocity(particle::Array{String,1}, global_best::Array{String,1}, local_best::Array{String,1}, alpha::Float64, beta::Float64)
    #@printf("current_position: %s\n", string(particle))
    #@printf("local_best: %s\n", string(local_best))
    #@printf("global_best: %s\n", string(global_best))
    swaps_to_local_best = get_swap_sequence(particle, local_best)
    swaps_to_global_best = get_swap_sequence(particle, global_best)
    @printf("swaps_to_global_best: %s\n", swaps_to_global_best)
    @printf("swaps_to_local_best: %s\n", swaps_to_local_best)
    velocity = []

    for lswap in swaps_to_local_best
        r = Random.rand()
        if r < alpha
            push!(velocity, lswap)
        end
    end
    for gswap in swaps_to_global_best
        r = Random.rand()
        if r < beta
            push!(velocity, gswap)
        end
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

    aligned_sequences = Array{String,1}(undef,length(sequences))
    predecessors = Array{Array{Int64,1}}(undef,0) # store which sequence indexes are at or below each other sequence in the tree

    sets = Array{Array{Int64,1},1}(undef,0)

    for i in 1:length(sequences)
        #seq in sequences
        seq = sequences[i]
        aligned_sequences[i] = seq
        #push!(aligned_sequences, deepcopy(seq))
        push!(predecessors, [])
        push!(sets, [i])
    end
    #@printf("Sequences: %s\n", string(sequences))
    @printf("Edges: %s\n", string(edges))


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
        dprintf("\n\nPerforming global alignment of %d %d\n%s\nand\n%s\nedge_weight: %d\n", idxA, idxB, A, B, min_edge[3])
        score, alignedA, alignedB = global_align(A, B)
        dprintf("Got alignment\n%s\nand\n%s\n", alignedA, alignedB)

        # TODO take new_gap_indexes_A, and insert same gaps into all current
        # predecessors of A.  And same for B
        new_gap_indexes_A = get_new_gap_indexes(A, alignedA)
        #@printf("new_gap_indexes_A: %s\n", string(new_gap_indexes_A))
        #@printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))

        #seqs_to_update = copy(sets[idxA])
        #deleteat!(seqs_to_update, findfirst(x -> x == idxA, seqs_to_update))
        seqs_to_update = [x for x in sets[idxA] if x != idxA]

        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_A)
        #@printf("predecessors A(%d): %s\n", idxA, string(sets[idxA]))

        new_gap_indexes_B = get_new_gap_indexes(B, alignedB)
        #@printf("new_gap_indexes_B: %s\n", string(new_gap_indexes_B))
        #println()
        #@printf("predecessors B(%d): %s\n", idxB, string(sets[idxB]))

        seqs_to_update = [x for x in sets[idxA] if x != idxA]
        #seqs_to_update = copy(sets[idxB])
        #deleteat!(seqs_to_update, findfirst(x -> x == idxB, seqs_to_update))

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
    sequences::Array{String,1}              # initial order/position
    aligned_sequences::Array{String,1}      # hold the aligned sequences to avoid unnecessary recomputes
    #sequence_map::Array{SequenceMap,1}
    position::Array{String,1} # ordering    # current order/position
    best_sequences::Array{String,1}         # best order/position
    best_score::Int                         # best score
    pidx::Int                               # particle index (for reverse mapping into the particle array)
    edges::Array{Tuple{String,String,Int},1}# edges for the current position
end

function print_particle(p::Particle)
    print("Particle: ")
    println(p.pidx)
    print("  position:\n")
    for seq in p.position
        print("    ")
        println(seq)
    end
    print("  best_sequences:\n")
    for seq in p.best_sequences
        print("    ")
        println(seq)
    end
    print("  best_score: ")
    println(p.best_score)
    print("  edges:")
    for edge in p.edges
        print(" ")
        print(edge[1])
    end
    println()
end


function PSO_MSA(sequences::Array{String,1}, iterations::Int)
    # N = 40 # length of the generated sequences
    # t = 10 # number of sequences
    # iterations = 2000
    t = length(sequences)
    solution_space = factorial(big(t)) #/ (t*(t-1))

    num_particles = 0
    if t > 10
        num_particles = 100
    else
        num_particles = t*(t-1) # num edges
    end
    @printf("num_particles: %d\n", num_particles)
    search_space = num_particles * iterations
    search_solution_ratio = search_space / solution_space

    @printf("Solution Space: %.02f\n", solution_space)
    @printf("Search Space: %.02f\n", search_space)
    @printf("Search/Solution space: %.02f\n", search_solution_ratio)

    nodes = Array{String,1}(undef,0)
    for i in 1:length(sequences)
        push!(nodes, sequences[i])
    end

    permutations = random_permutations(nodes, num_particles)

    # Initialize particles list from the generated permutations
    println("Initializing particles")
    particles = Array{Particle,1}(undef,length(permutations))
    for pidx in 1:length(permutations)
        sequences = copy(permutations[pidx])
        aligned_sequences = copy(sequences)

        p_edges = Array{Tuple{String,String,Int},1}(undef,0)
        for sidx in 2:length(sequences)
            seqA = sequences[sidx-1]
            seqB = sequences[sidx]
            scoreAB = score_sequences([seqA, seqB])
            push!(p_edges, (seqA, seqB, scoreAB))
        end

        aligned_sequences = progressive_alignment_inorder(sequences, p_edges)
        initial_best_score = score_sequences(aligned_sequences)
        p = Particle(
            permutations[pidx],         # original sequences/position
            aligned_sequences,          # aligned_sequences
            copy(permutations[pidx]),   # current_position
            copy(permutations[pidx]),   # best_sequences
            initial_best_score,         # initial best score
            pidx,                       # index
            p_edges)                    # edges for the current position
        particles[pidx] = p
        print_particle(particles[pidx])
    end

    #alpha = Random.rand()
    #beta = Random.rand()
    alpha = 0.2
    beta = 0.2

    # keep track of the local best score+ordering for each particle
    local_best_scores = [p.best_score for p in particles]
    #local_best_particles = deepcopy(particles)
    # keep track of the global best particle
    global_best_particle = particles[argmax(local_best_scores)]
    println("Initial Best")
    @printf("Global best (score = %d, pidx = %d):\n", global_best_particle.best_score, global_best_particle.pidx)
    for seq in global_best_particle.sequences
        println(seq)
    end
    @printf("Initial global best score: %d\n", global_best_particle.best_score)

    # Pick a starting position, arbitrary
    #position_xid = copy(particles[1].sequences)

    iteration_metrics = []

    random_subgroup_move_max = 10
    random_subgroup_move_count = 0

    # percentage of particles which didn't move
    reshuffle_threshold = length(particles) * 0.75

    for iteration in 1:iterations
        # Revert all particle aligned sequences to their original unaligned values
        # for p in particles
        #     p.aligned_sequences = copy(p.position)
        # end

        unmoved_particles = []

        @printf("Beginning iteration %d\n", iteration)
        start = Dates.now()

        if random_subgroup_move_count >= random_subgroup_move_max
            do_random_subgroup_move = true
            random_subgroup_move_count = 0
        else
            random_subgroup_move_count += 1
            do_random_subgroup_move = false
        end

        for pidx in 1:length(particles)
            particle = particles[pidx]
            println("Looking at particle:")
            print_particle(particles[pidx])

            # dprintf("Global best (score=%d):\n", global_best_particle.best_score)
            # for seq in global_best_particle.position
            #     dprintf("  ")
            #     dprintf("%s\n", seq)
            # end
            new_velocity = nothing

            if do_random_subgroup_move
                random_subgroup_prob = 0.1
                if Random.rand() < random_subgroup_prob
                    random_position = random_permutations(particle.position, 1)[1]
                    @printf("Randomly moving particle from %s to %s\n", string(particle.position), string(random_position))
                    new_velocity = pso_particle_velocity(
                        random_position,
                        global_best_particle.best_sequences,
                        particle.best_sequences,
                        alpha, beta)
                end
            end
            if new_velocity === nothing
                new_velocity = pso_particle_velocity(
                    particle.position,
                    global_best_particle.best_sequences,
                    particle.best_sequences,
                    alpha, beta)
            end


            if length(new_velocity) == 0
                push!(unmoved_particles, pidx)
            end

            if length(new_velocity) > 0
                # Apply the velocity to the position of the particle
                new_position = apply_velocity_to_particle(particle.position, new_velocity)
                @printf("Applied velocity %s to particle %d, moved from\n%s, to %s\n",
                    string(new_velocity),
                    particle.pidx,
                    particle.position,
                    new_position)
                particle.position = new_position

                # Calculate the new edges for the current position (sequence ordering)
                p_edges = Array{Tuple{String,String,Int},1}(undef,length(particle.position)-1)
                for sidx in 2:length(particle.position)
                    seqA = particle.position[sidx-1]
                    seqB = particle.position[sidx]
                    scoreAB = score_sequences([seqA, seqB])
                    p_edges[sidx-1] = (seqA, seqB, scoreAB)
                    #push!(p_edges, (seqA, seqB, scoreAB))
                end
                particle.edges = p_edges

                # Perform a progressive alignment of the new position with the new edges
                aligned_sequences = progressive_alignment_inorder(particle.position, particle.edges)
                score = score_sequences(aligned_sequences)

                particle.aligned_sequences = copy(aligned_sequences)

                # if the new position score is better, update the local best
                if score > particle.best_score
                    @printf("Updating particle best score to %d\n", score)
                    particle.best_score = score
                    particle.best_sequences = copy(particle.position)
                end

                # if particle.best_score > global_best_particle.best_score
                #     @printf("Updating global best score to %d\n", particle.best_score)
                #     global_best_particle = deepcopy(particle)
                # end
            end

            # Update the global best at the end of the round
            particle_best_scores = [p.best_score for p in particles]
            best_particle = particles[argmax(particle_best_scores)]
            if best_particle.best_score > global_best_particle.best_score
                @printf("Updating global best score to %d\n", best_particle.best_score)
                global_best_particle = deepcopy(best_particle)
            end
        end

        if length(unmoved_particles) >= 0
            @printf("Reshuffling particles: %s\n", string(unmoved_particles))
            for reshuffle_pidx in unmoved_particles
                p = particles[reshuffle_pidx]
                p.position = random_permutations(p.position, 1)[1]
            end
        end

        end_time = Dates.now()
        duration = end_time - start
        push!(iteration_metrics, duration)
        @printf("\nIteration %d took %d milliseconds\n", iteration, duration.value)
        @printf("Global best (score = %d, pidx = %d):\n",
            global_best_particle.best_score,
            global_best_particle.pidx)
        for seq in global_best_particle.aligned_sequences
            println(seq)
        end
        println()
    end

    println("FINISHED PARTICLE SWARM")

    # print the final results
    @printf("\nGlobal best (score = %d, pidx = %d):\n",
        global_best_particle.best_score,
        global_best_particle.pidx)
    for seq in global_best_particle.aligned_sequences
        println(seq)
    end
    println()

    # # remove any all-gap columns
    # i = 1
    # while i <= length(global_best_particle.aligned_sequences[1])
    # #for i in 1:length(global_best_particle_value[1])
    #     #@printf("i = %d\n", i)
    #     all_gaps = true
    #     for seq in global_best_particle.aligned_sequences
    #         if seq[i] != '-'
    #             all_gaps = false
    #             break
    #         end
    #     end
    #     if all_gaps
    #         println("removing gap-only column")
    #         # remove the column at index i
    #         for seqidx in 1:length(global_best_particle.aligned_sequences)
    #             seq = global_best_particle.aligned_sequences[seqidx]
    #             seq = string(seq[1:i-1], seq[i+1:end])
    #             global_best_particle.aligned_sequences[seqidx] = seq
    #         end
    #         i -= 1
    #         #println(i)
    #     end
    #     i += 1
    # end
    # global_best_particle.best_score = score_sequences(global_best_particle.aligned_sequences)
    # @printf("\nGlobal best AFTER REMOVING GAP COLUMNS (score = %d):\n", global_best_particle.best_score)
    # for seq in global_best_particle.aligned_sequences
    #     println(seq)
    # end
    # println()
end

function main()
    t = 5 # number of sequences
    N = 15
    sequences = generate_sequences(t, N)
    println(sequences)
    perf = @timed PSO_MSA(sequences, 200)

    println(perf)

    # swap_sequence = get_swap_sequence([1,2,3,4,5,6], [4,3,2,1,6,5])
    # println(swap_sequence)
end
main()