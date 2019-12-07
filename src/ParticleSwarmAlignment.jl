module ParticleSwarmAlignment

export PSO_MSA, score_sequences, generate_sequences, global_align, progressive_alignment_inorder

using Random
using Printf
using Dates

# Random module must be seeded with a constant in order to have
# reproducible test results
Random.seed!(1)

# Debug printf function which can be turned on/off with debug flag
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

# Needleman-Wunsch Global Alignment using fixed penalties.
# Extending with keyword arguments for optional preallocated s, b matrices
# will drastically improve performance during its use in progressive alignment
function global_align(v, w, match_penalty=1, mismatch_penalty=-1, deletion_penalty=-1)
    n1 = length(v)
    n2 = length(w)
    #if !use_preallocated_matrices
    s = zeros(Float64, n1+1, n2+1)
    b = zeros(Float64, n1+1, n2+1)
    #end

    for i in 1:(n1+1)
        s[i,1] = (i-1) * deletion_penalty
        b[i,1] = 2
    end
    for j in 1:(n2+1)
        s[1,j] = (j-1) * deletion_penalty
        b[1,j] = 3
    end

    for i in 2:(n1+1)
        for j in 2:(n2+1)
            if v[i-1] == w[j-1]
                ms = s[i-1,j-1] + match_penalty
            else
                # ignore cases where a letter is paired with a gap
                # do not consider this a mismatch
                # if v[i-1] != '-' && w[j-1] != '-'
                #     ms = s[i-1,j-1] + mismatch_penalty
                # else
                #     # if a letter is paired with a gap, add no penalty
                #     ms = s[i-1,j-1] #+ match_penalty # + 0.5 * mismatch_penalty
                # end
                ms = s[i-1,j-1] + mismatch_penalty
            end
            test = [ms, s[i-1,j] + deletion_penalty, s[i,j-1] + deletion_penalty]
            p = argmax(test)
            s[i,j] = test[p]
            b[i,j] = p
        end
    end

    i = n1+1
    j = n2+1
    sv = []
    sw = []
    while(i > 1 || j > 1)
        p = b[i,j]
        if (p == 1)
            i = i-1
            j = j-1
            push!(sv, v[i])
            push!(sw, w[j])
        elseif p == 2
            i=i-1
            push!(sv, v[i])
            push!(sw, "-")
        elseif p == 3
            j = j-1
            push!(sv, "-")
            push!(sw, w[j])
        else
            break
        end
    end

    return (s[n1+1,n2+1], join(reverse(sv)), join(reverse(sw)))
end


# Generates t sequences all of length l, and returns an array of strings
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


# Performs an *ordered* equals comparison of two ordered arrays
function array_equals_ordered(A::Array, B::Array)
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

# Determines whether Any element e is contained within Any array A.
# Returns true if so.
function array_contains(A::Array, e)
    for i in 1:length(A)
        if array_equals_ordered(A[i], e)
            return true
        end
    end
    return false
end


# Creates num random permutations of the elements of A. If must_be_unique is left
# as true, then it ensures that the permutations returned are unique. This can fail
# if it takes too long to find a unique permutation, as the method randomly searches
# for permutations instead of iterating through them in order.
function random_permutations(A::Array{String,1}, num::Int64, must_be_unique::Bool=true)::Array{Array{String,1},1}
    # Return num random permutations of the elements of A
    results = Array{Array{String,1},1}(undef,0)
    if num >= factorial(big(length(A)))
        throw(ErrorException("num is too high"))
    end
    max_check = 100
    try_count = 0
    for i in 1:num
        localA = deepcopy(A)
        randomA = shuffle!(localA)
        if must_be_unique
            while array_contains(results, randomA)
                randomA = shuffle!(A)
                try_count += 1
                if try_count > max_check
                    msg = @sprintf("Exceeded max tries (%d) while trying to find unique permutations\n", max_check)
                    throw(ErrorException(msg))
                end
            end
        end
        push!(results, randomA)
    end
    return results
end

# Returns a swap sequence [(idx1, idx2), ...] to turn A into B
# Not necessarily the optimal swap sequence, of the least number of swaps
function get_swap_sequence(A_in::Array{String,1}, B_in::Array{String,1})
    if length(A_in) != length(B_in)
        throw(ErrorException("Length of A and B must be the same"))
    end
    dprintf("Getting swaps between:\n%s,\n%s\n", string(A_in), string(B_in))
    A = Array{String,1}(undef,length(A_in))
    B = Array{String,1}(undef,length(B_in))
    for i in 1:length(A_in)
        A[i] = A_in[i]
        B[i] = B_in[i]
    end

    swap_sequence = Tuple{Int,Int}[]

    for i in 1:length(A)
        if A[i] != B[i]
            # We want to turn A into B
            # so find where B[i] is in A[i:end], and swap that element in A with the element at A[i]
            a_index = findfirst(x -> x == B[i], A[i:end])
            if a_index === nothing
                msg = string("Could not find ", string(B[i], string(" in ", string(A[i:end]))))
                throw(ErrorException(msg))
            end
            a_index += i - 1

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


# Calculates a new "velocity" for the particle position given global best position,
# and the local best position. Alpha and Beta are tuning parameters which determine
# how quickly the particle trends towards the local/global best.
function pso_particle_velocity(particle::Array{String,1}, global_best::Array{String,1}, local_best::Array{String,1}, alpha::Float64, beta::Float64)
    dprintf("current_position: %s\n", string(particle))
    dprintf("local_best: %s\n", string(local_best))
    dprintf("global_best: %s\n", string(global_best))
    swaps_to_local_best = get_swap_sequence(particle, local_best)
    swaps_to_global_best = get_swap_sequence(particle, global_best)
    dprintf("swaps_to_global_best: %s\n", swaps_to_global_best)
    dprintf("swaps_to_local_best: %s\n", swaps_to_local_best)
    velocity = Tuple{Int,Int}[]

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


# This takes a particle position, and performs the swaps specified by
# swap_sequence_velocity, in order. Returns the resulting position array.
function apply_velocity_to_particle(position_in::Array{String,1}, swap_sequence_velocity::Array)::Array{String,1}
    position = deepcopy(position_in)
    for (idxA, idxB) in swap_sequence_velocity
        temp = position[idxA]
        position[idxA] = position[idxB]
        position[idxB] = temp
    end
    return position
end


# Determines where gaps need to be inserted into old to obtain the new string value.
# Returns the indexes.
function get_new_gap_indexes(old::String, new::String)::Array{Int64,1}
    old_local = old
    indexes = Int[]
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

# Determines where gaps need to be inserted into old to obtain the new string value.
# Returns the indexes.
# This is an optimization of get_new_gap_indexes which does incur large memory copies
# due to string copying.
function get_new_gap_indexes2(old::String, new::String)::Array{Int64,1}
    indexes = Int[]
    old_i = 1
    new_i = 1
    while new_i <= length(new)
        # if new string had a gap either matched with an old non-gap, or past the end of the old string
        if (new[new_i] == '-' && (old_i > length(old) || old[old_i] != '-'))
            # new gap
            push!(indexes, new_i)
            new_i += 1
        elseif (new[new_i] == old[old_i])
            new_i += 1
            old_i += 1
        else
            @printf("Unknown error old[%d] = %s, new[%d] = %s\n", old_i, string(old[old_i]), new_i, string(new[new_i]))
            throw(ErrorException("unknown error!"))
        end
    end
    return indexes
end

# Inserts gap characters into the sequences specified by the indexes in seq_indexes.
# The gaps are inserted at the indexes within each specified sequence at positions gap_indexes.
function insert_gaps_at(sequences::Array, seq_indexes::Array{Int64,1}, gap_indexes::Array{Int64,1})
    if length(seq_indexes) > 0 && length(gap_indexes) > 0
        gap_indexes = sort(gap_indexes)

        for sidx in 1:length(seq_indexes)
            seq_index = seq_indexes[sidx]
            s = sequences[seq_index]

            for i in gap_indexes
                # insert a '-' character at index i
                s = string(string(s[1:i-1], '-'), s[i:end])
            end
            # put the modified string back in the array
            sequences[seq_index] = s
        end
    end
end

# k is the length of the sequences, t is the number of sequences, and A is the sequences.
# The sequences of A must be the same length, or k must be the length of the longest sequence.
function Make_Profile(k::Int, t::Int, A::Array{String,1})::Dict{String,Array{Int64,1}}
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

# Score the sequences in A by obtaining the consensus score using Make_Profile
function score_sequences(A::Array{String,1})::Int
    k = length(A[1])    # how long are the sequences
    t = length(A)       # how many sequences
    profile = Make_Profile(k, t, A)
    # score is the sum of the maximum occurring letter (not including gaps) in each position
    score = 0
    for i in 1:k
        a = profile["A"][i]
        t = profile["T"][i]
        g = profile["G"][i]
        c = profile["C"][i]
        score += max(a, t, g, c)
    end
    return score
end

# Performs a comparison of A and B, but ignoring the gap characters in each.
function sequence_equals_ignore_gaps(A::String, B::String)
    A = filter(x -> x != '-', A)
    B = filter(x -> x != '-', B)
    return A == B
end


# This performs a progressive alignment of sequences, using the elements of the edges array
# to determine the order with which to perform neighbor alignments.
# The edges array must contain the pairwise edges and weights for all consecutive
# elements of the sequences array.
function progressive_alignment_inorder(sequences::Array, edges::Array{Tuple{String,String,Int},1})::Array{String,1}
    if length(edges) != length(sequences) - 1
        msg = @sprintf("Expected %d edge weights, got %d", length(sequences) - 1, length(edges))
        throw(ErrorException(msg))
    end

    aligned_sequences = Array{String,1}(undef,length(sequences))

    sets = Array{Array{Int64,1},1}(undef,0)

    for i in 1:length(sequences)
        #seq in sequences
        seq = sequences[i]
        aligned_sequences[i] = seq
        push!(sets, [i])
    end
    dprintf("Edges: %s\n", string(edges))

    local_edges = sort(edges, by=x -> x[3], rev=true) # sort desc by weight

    for i in 1:length(edges)
        min_edge = pop!(local_edges)

        dprintf("min_edge: %s\n", string(min_edge))

        # use aligned sequences here instead of original sequences
        #idxA = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[1]), sequences)
        idxA = findfirst(s -> s == min_edge[1], sequences)
        #idxB = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[2]), sequences)
        idxB = findfirst(s -> s == min_edge[2], sequences)
        A = aligned_sequences[idxA]
        B = aligned_sequences[idxB]
        dprintf("\n\nPerforming global alignment of %d %d\n%s\nand\n%s\nedge_weight: %d\n", idxA, idxB, A, B, min_edge[3])
        score, alignedA, alignedB = global_align(A, B)
        dprintf("Got alignment\n%s\nand\n%s\n", alignedA, alignedB)

        # take new_gap_indexes_A, and insert same gaps into all current
        # predecessors of A.  And same for B
        new_gap_indexes_A = get_new_gap_indexes2(A, alignedA)
        seqs_to_update = [x for x in sets[idxA] if x != idxA]
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_A)

        new_gap_indexes_B = get_new_gap_indexes2(B, alignedB)
        seqs_to_update = [x for x in sets[idxB] if x != idxB]
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_B)

        aligned_sequences[idxA] = alignedA
        aligned_sequences[idxB] = alignedB

        # Update predecessors for A and B
        # TODO avoid the loop somehow by reusing the sets array from one of the existing idxA or idxB elements
        newset = Int[]#[idxA, idxB]
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


# Searches the edges array to find the edge between nodeA and nodeB.
# Returns the edge Tuple{String,String,Int}, as well as the index it was found at.
function get_edge_between(nodeA::String, nodeB::String, edges::Array)
    #nodeA = filter(x -> x != '-', nodeA)
    #nodeB = filter(x -> x != '-', nodeB)
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

# Holds all relevant information for a particle in the swarm
mutable struct Particle
    sequences::Array{String,1}              # initial order/position
    aligned_sequences::Array{String,1}      # hold the aligned sequences to avoid unnecessary recomputes
    position::Array{String,1} # ordering    # current order/position
    best_sequences::Array{String,1}         # best order/position
    best_score::Int                         # best score
    pidx::Int                               # particle index (for reverse mapping into the particle array)
    edges::Array{Tuple{String,String,Int},1}# edges for the current position
end

# Utility function for printing the more relevant information about a particle for debugging.
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


# The primary function of the module, which takes the sequences, the number of swarm iterations to run, and
# optionally the number of particles to use. If no num_particles is provided, then a reasonable value is
# calculated based on the number of sequences.
# The length of the sequences do not all need to be the same, however the performance may reduce if they
# are not because additional work must be done to initially align neighbors before computing edge weights
function PSO_MSA(sequences::Array{String,1}, iterations::Int, num_particles::Int=0)
    # Number of sequences
    t = length(sequences)
    # Number of permutations ignoring complete reversals
    solution_space = factorial(big(t)) / 2

    if num_particles == 0
        if t >= 8
            num_particles = 100
        else
            #convert(Int64, t*(t-1)/2) # num edges over 2
            num_particles = t
        end
    end
    @printf("Number of Particles: %d\n", num_particles)
    search_space = num_particles * iterations
    search_solution_ratio = search_space / solution_space

    @printf("Solution Space: %.02f\n", solution_space)
    @printf("Search Space: %.02f\n", search_space)
    @printf("Search/Solution space: %.010f\n", search_solution_ratio)

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
        @printf("Initializing particle %d\n", pidx)
        # Initial edges
        p_edges = Array{Tuple{String,String,Int},1}(undef,0)
        for sidx in 2:length(sequences)
            seqA_orig = sequences[sidx-1]
            seqB_orig = sequences[sidx]
            #scoreAB = score_sequences([seqA, seqB])
            # if length(seqA_orig) != length(seqB_orig)
            #     (align_score, seqA, seqB), t, bytes, gct, memalloc = @timed global_align(seqA_orig, seqB_orig)
            #     println(t)
            # else
            #     seqA, seqB = seqA_orig, seqB_orig
            # end
            # scoreAB = score_sequences([seqA, seqB])

            # Here the initial score for each neighbor distance is set to zero
            # to avoid having to align them to obtain a valid distance score
            scoreAB = 0
            push!(p_edges, (seqA_orig, seqB_orig, scoreAB))
        end

        @time aligned_sequences = progressive_alignment_inorder(sequences, p_edges)

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
        if debug
            print_particle(particles[pidx])
        end
    end

    #alpha = Random.rand()
    #beta = Random.rand()

    # Here are two important tuning values. These determine how rapidly the particles trend
    # towards their own past best and the swarm's best solution. With some non-systematic testing
    # it was determined that 0.2 is a reasonable value for both. Too high of a value results
    # in particles not searching the area around the local optimum enough. Too low of a value results
    # in not enough motion happening.
    alpha = 0.2
    beta = 0.2

    # keep track of the local best score+ordering for each particle
    local_best_scores = [p.best_score for p in particles]
    # keep track of the global best particle
    global_best_particle = particles[argmax(local_best_scores)]
    println("Initial Best")
    @printf("Global best (score = %d, pidx = %d):\n", global_best_particle.best_score, global_best_particle.pidx)
    for seq in global_best_particle.sequences
        println(seq)
    end
    @printf("Initial global best score: %d\n", global_best_particle.best_score)

    iteration_metrics = Millisecond[]

    # These take a potential solution to swarm stagnation. Every so often, we can reassign
    # a certain number of particles in the swarm to a random position.
    # Setting random_subgroup_prob to -1 disables this. We have disabled it due to
    # it introducing too much randomness into the system.
    random_subgroup_move_max = 10
    random_subgroup_move_count = 0
    random_subgroup_prob = -1

    # Percentage of particles which if they don't move, move them randomly
    # This is also a significant tuning value. Setting it to one means that any time a particle
    # doesn't move for the past iteration, immediately pick it up and give it a new position.
    #reshuffle_threshold = length(particles) * 0.5
    reshuffle_threshold = 1

    for iteration in 1:iterations

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
            if debug
                println("Looking at particle:")
                print_particle(particles[pidx])
            end

            new_velocity = nothing

            if do_random_subgroup_move
                if Random.rand() < random_subgroup_prob
                    random_position = random_permutations(particle.position, 1)[1]
                    dprintf("Randomly moving particle from %s to %s\n", string(particle.position), string(random_position))
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

            # This particle did not move in this iteration
            if length(new_velocity) == 0
                push!(unmoved_particles, pidx)
            end

            if length(new_velocity) > 0
                # Apply the velocity to the position of the particle
                new_position = apply_velocity_to_particle(particle.position, new_velocity)
                dprintf("Applied velocity %s to particle %d, moved from\n%s, to %s\n",
                    string(new_velocity),
                    particle.pidx,
                    particle.position,
                    new_position)
                particle.position = new_position

                # Perform the same swaps to the aligned sequences
                particle.aligned_sequences = apply_velocity_to_particle(particle.aligned_sequences, new_velocity)

                # Calculate the new edges for the current position (sequence ordering)
                p_edges = Array{Tuple{String,String,Int},1}(undef,length(particle.position)-1)
                for sidx in 2:length(particle.position)
                    seqA_orig = particle.position[sidx-1]
                    seqB_orig = particle.position[sidx]
                    seqA_aligned = particle.aligned_sequences[sidx-1]
                    seqB_aligned = particle.aligned_sequences[sidx]

                    scoreAB = score_sequences([seqA_aligned, seqB_aligned])
                    p_edges[sidx-1] = (seqA_orig, seqB_orig, scoreAB)
                end
                particle.edges = p_edges

                # Perform a progressive alignment of the new position with the new edges
                aligned_sequences = progressive_alignment_inorder(particle.position, particle.edges)
                score = score_sequences(aligned_sequences)

                particle.aligned_sequences = aligned_sequences

                # if the new position score is better, update the local best
                if score > particle.best_score
                    dprintf("Updating particle best score to %d\n", score)
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

        if length(unmoved_particles) >= reshuffle_threshold
            dprintf("Reshuffling particles: %s\n", string(unmoved_particles))
            for reshuffle_pidx in unmoved_particles
                p = particles[reshuffle_pidx]
                p.position = random_permutations(p.position, 1)[1]
            end
        end

        end_time = Dates.now()
        duration = end_time - start
        push!(iteration_metrics, duration)
        @printf("Iteration %d took %d milliseconds\n", iteration, duration.value)
        if debug
            @printf("Global best (score = %d, pidx = %d):\n",
                global_best_particle.best_score,
                global_best_particle.pidx)
            for seq in global_best_particle.aligned_sequences
                println(seq)
            end
            println()
        end
    end

    println("FINISHED PARTICLE SWARM")

    # print the final results
    @printf("Global best (score = %d, pidx = %d):\n",
        global_best_particle.best_score,
        global_best_particle.pidx)
    for seq in global_best_particle.aligned_sequences
        println(seq)
    end
    println()

    return global_best_particle.best_score, global_best_particle.aligned_sequences
end

end # end module