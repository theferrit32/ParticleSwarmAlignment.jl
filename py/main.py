import random
import time
import math
import copy
import datetime

# local import
from align import global_align

random.seed(1)
debug = False

def dprintf(s, *args):
    if debug:
        print(s % (args), end='')

def Make_Profile(k:int, t:int, A:list) -> dict:
    idxA = [0] * k
    idxC = [0] * k
    idxG = [0] * k
    idxT = [0] * k
    idxGap = [0] * k

    for index in range(0, k):
        for seq in A:
            if seq[index] == 'A':
                idxA[index] += 1
            elif seq[index] == 'C':
                idxC[index] += 1
            elif seq[index] == 'G':
                idxG[index] += 1
            elif seq[index] == 'T':
                idxT[index] += 1
            elif seq[index] == '-':
                idxGap[index] += 1
            else:
                msg = "character [%s] not allowed in DNA sequence" % str(seq[index])
                raise RuntimeError(msg)

    return {
        "A": idxA,
        "C": idxC,
        "G": idxG,
        "T": idxT,
        "-": idxGap
    }

def score_sequences(A:list) -> int:
    k = len(A[1])    # how long are the sequences
    t = len(A)       # how many sequences
    profile = Make_Profile(k, t, A)
    # score is the sum of the maximum occurring letter (not including gaps) in each position
    score = 0
    #for i in 1:k
    for i in range(0, k):
        score += max(profile["A"][i], profile["T"][i], profile["G"][i], profile["C"][i])
    return score

def generate_sequences(t:int, l:int) -> list:
    # t is the number of sequences to create
    # l is the length of the sequences
    DNA = []
    base_arr = ["A", "T", "G", "C"]

    for t_index in range(0, t):
        DNA.append("")
        for _ in range(0, l):
            r = int(math.floor(random.random() * 4))
            DNA[t_index] = DNA[t_index] + base_arr[r]

    return DNA


def array_equals_ordered(A:list, B:list) -> bool:
    #dprintf("Comparing %s to %s\n", string(A), string(B))
    if len(A) != len(B):
        return False

    for i in range(0, len(A)):
        if A[i] != B[i]:
            return False
    return True

def array_contains(A:list, e) -> bool:
    for i in range(0, len(A)):
        if array_equals_ordered(A[i], e):
            return True
    return False

def random_permutations(A:list, num:int) -> list:
    # Return num random permutations of the elements of A
    # results = Array{Array{String,1},1}(undef,0)
    # if num >= factorial(big(length(A)))
    #     throw(ErrorException("num is too high"))
    # end
    # max_check = 50
    # try_count = 0
    # for i in 1:num
    #     localA = deepcopy(A)
    #     randomA = shuffle!(localA)
    #     while array_contains(results, randomA)
    #         randomA = shuffle!(A)
    #         try_count += 1
    #         if try_count > max_check
    #             msg = @sprintf("Exceeded max permutation tries: %d", max_check)
    #             throw(ErrorException(msg))
    #         end
    #     end
    #     push!(results, randomA)
    # end
    # return results
    results = []
    if num >= math.factorial(len(A)):
        raise RuntimeError("num is too high: %s, max is: %s" % (num, math.factorial(len(A))))
    max_try = 50
    try_count = 0
    for _ in range(0, num):
        localA = copy.deepcopy(A)
        random.shuffle(localA)
        while array_contains(results, localA):
            random.shuffle(localA)
            try_count += 1
            if try_count > max_try:
                raise RuntimeError("Exceeded max permutation tries: %d" % max_try)
        results.append(localA)
    #print("random_permutations: %s" % str(results))
    return results


# Returns a swap sequence [(idx1, idx2), ...] to turn A into B
# Not necessarily the optimal swap sequence
def get_swap_sequence(A_in:list, B_in:list) -> list:
    if len(A_in) != len(B_in):
        raise RuntimeError("Length of A and B must be the same")

    #@printf("Getting swaps between:\n%s,\n%s\n", string(A_in), string(B_in))
    #A = Array{String,1}(undef,length(A_in))
    #B = Array{String,1}(undef,length(B_in))
    # for i in 1:length(A_in)
    #     A[i] = A_in[i]
    #     B[i] = B_in[i]
    # end
    A = []
    B = []
    for _ in range(len(A_in)):
        A.append(None)
        B.append(None)

    for i in range(0, len(A_in)):
        A[i] = A_in[i]
        B[i] = B_in[i]

    swap_sequence = []

    #for i in 1:length(A)
    for i in range(0, len(A)):
        #@printf("\nChecking index %d,\nlocalA: %s\nlocalB: %s\n", i, string(A), string(B))
        dprintf("Checking index %d\nA = %s\nB = %s" % (i, str(A), str(B)))
        if A[i] != B[i]:
            # We want to turn A into B
            # so find where B[i] is in A[i:end], and swap that element in A with the element at A[i]
            #a_index = findfirst(x -> x == B[i], A[i:end])

            if B[i] not in A[i:]:
                raise RuntimeError("Could not find %s in %s" % (B[i], A[i:]))

            a_index = A[i:].index(B[i])
            a_index += i

            dprintf("Doing swap in A (%d, %d)\n", i, a_index)
            temp = A[i]
            A[i] = A[a_index]
            A[a_index] = temp

            #push!(swap_sequence, (i, a_index))
            swap_sequence.append((i, a_index))

    equals_check = array_equals_ordered(A, B)
    if equals_check is False:
        raise RuntimeError("A and B were not the same:\n%s\n%s\n", str(A), str(B))

    return swap_sequence

# Calculates a new "velocity" for the particle given global best particle, local best particle
def pso_particle_velocity(particle:list, global_best:list, local_best:list, alpha:float, beta:float) -> list:
    swaps_to_local_best = get_swap_sequence(particle, local_best)
    swaps_to_global_best = get_swap_sequence(particle, global_best)
    #dprintf("swaps_to_global_best: %s\n", swaps_to_global_best)
    #dprintf("swaps_to_local_best: %s\n", swaps_to_local_best)
    velocity = []

    for lswap in swaps_to_local_best:
        r = random.random()
        if r < alpha:
            velocity.append(lswap)

    for gswap in swaps_to_global_best:
        r = random.random()
        if r < beta:
            velocity.append(gswap)

    return velocity


def apply_velocity_to_particle(particle_in:list, swap_sequence_velocity:list) -> list:
    particle = copy.deepcopy(particle_in)
    for (idxA, idxB) in swap_sequence_velocity:
        temp = particle[idxA]
        particle[idxA] = particle[idxB]
        particle[idxB] = temp
    return particle

def get_new_gap_indexes(old:str, new:str) -> list:
    old_local = old
    #old_array = split(old, "")
    indexes = []
    dprintf("get_new_gap_indexes, old = %s, new = %s" % (old, new))
    for i in range(0, len(new)):
        if (new[i] == '-' and i == (len(old_local))) or (new[i] == '-' and old_local[i] != '-'):
            #@printf("get_new_gap_indexes, new gap at %d, old: %s, new: %s\n", i, old_local, new)
            #push!(indexes, i)
            #old_local = string(string(old_local[1:i-1], '-'), old_local[i:end])
            indexes.append(i)
            old_local = old_local[0:i] + "-" + old_local[i:]
        elif new[i] == old_local[i]:
            continue
        else:
            raise RuntimeError("Sequences ignoring gaps were not the same")

    return indexes


def get_new_gap_indexes2(old:str, new:str) -> list:
    indexes = []
    old_i = 0
    new_i = 0
    while new_i <= len(new):
        # if new string had a gap either matched with an old non-gap, or past the end of the old string
        if (new[new_i] == '-' and (old_i > len(old) or old[old_i] != '-')):
            # new gap
            #push!(indexes, new_i)
            indexes.append(new_i)
            new_i += 1
        elif (new[new_i] == old[old_i]):
            new_i += 1
            old_i += 1
        else:
            #@printf("Unknown error old[%d] = %s, new[%d] = %s\n", old_i, string(old[old_i]), new_i, string(new[new_i]))
            #throw(ErrorException("unknown error!"))
            print("Unknown error old[%d] = %s, new[%d] = %s" % (old_i, old[old_i], new_i, new[new_i]))
            raise RuntimeError("unknown error!")

    return indexes


def insert_gaps_at(sequences:list, seq_indexes:list, gap_indexes:list):
    if len(seq_indexes) > 0 and len(gap_indexes) > 0:
        gap_indexes = sorted(gap_indexes)

        #for sidx in 1:length(seq_indexes)
        for sidx in range(0, len(seq_indexes)):
            seq_index = seq_indexes[sidx]
            s = sequences[seq_index]

            for i in gap_indexes:
                # insert a '-' character at index i
                #@printf("Inserting gap into %s at index %d\n", s, i)
                #s = string(string(s[1:i-1], '-'), s[i:end])
                s = s[0:i] + "-" + s[i:]
                #@printf("New value: %s\n", s)

            # put the modified string back in the array
            sequences[seq_index] = s


def sequence_equals_ignore_gaps(A:str, B:str):
    A = A.replace("-", "")
    B = B.replace("-", "")
    return A == B


def progressive_alignment_inorder(sequences:list, edges:list) -> list:
    if len(edges) != len(sequences) - 1:
        msg = "Expected %d edge weights, got %d" % (len(sequences) - 1, len(edges))
        raise RuntimeError(msg)

    aligned_sequences = [None] * len(sequences)
    sets = []

    for i in range(0, len(sequences)):
        #seq in sequences
        seq = sequences[i]
        aligned_sequences[i] = seq
        sets.append([i])

    dprintf("Edges: %s\n", str(edges))

    #local_edges = sort(edges, by=x -> x[3], rev=true) # sort desc by weight
    local_edges = sorted(edges, key=lambda x: x[2], reverse=True)

    for i in range(0, len(edges)):
        min_edge = local_edges.pop()

        #@printf("aligned_sequences: %s\n", string(aligned_sequences))
        #@printf("original edges: %s\n", string(edges))
        #@printf("local_edges: %s\n", string(local_edges))
        dprintf("min_edge: %s\n", str(min_edge))

        # use aligned sequences here instead of original sequences
        #idxA = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[1]), sequences)
        #idxA = findfirst(s -> s == min_edge[1], sequences)
        idxA = sequences.index(min_edge[0])
        #idxB = findfirst(s -> sequence_equals_ignore_gaps(s, min_edge[2]), sequences)
        #idxB = findfirst(s -> s == min_edge[2], sequences)
        idxB = sequences.index(min_edge[1])
        if idxA is None or idxB is None:
            raise RuntimeError("idxA = %s, idxB = %s" % (idxA, idxB))
        A = aligned_sequences[idxA]
        B = aligned_sequences[idxB]
        dprintf("\n\nPerforming global alignment of %d %d\n%s\nand\n%s\nedge_weight: %d\n", idxA, idxB, A, B, min_edge[2])
        score, alignedA, alignedB = global_align(A, B)
        dprintf("Got alignment\n%s\nand\n%s\n", alignedA, alignedB)

        # take new_gap_indexes_A, and insert same gaps into all current
        # predecessors of A.  And same for B
        new_gap_indexes_A = get_new_gap_indexes(A, alignedA)
        seqs_to_update = [x for x in sets[idxA] if x != idxA]
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_A)

        new_gap_indexes_B = get_new_gap_indexes(B, alignedB)
        seqs_to_update = [x for x in sets[idxB] if x != idxB]
        insert_gaps_at(aligned_sequences, seqs_to_update, new_gap_indexes_B)

        aligned_sequences[idxA] = alignedA
        aligned_sequences[idxB] = alignedB

        # for seq in aligned_sequences
        #     println(seq)
        # end

        # Update predecessors for A and B
        # TODO avoid the loop somehow by reusing the sets array from one of the existing idxA or idxB elements
        newset = []#[idxA, idxB]
        newset.extend(sets[idxA])
        newset.extend(sets[idxB])
        sets[idxA] = newset
        sets[idxB] = newset
        for p in newset:
            sets[p] = newset

    return aligned_sequences


# Return the (edge, index)
def get_edge_between(nodeA:str, nodeB:str, edges:list):
    #nodeA = filter(x -> x != '-', nodeA)
    #nodeB = filter(x -> x != '-', nodeB)

    for idx in range(0, len(edges)):
        edge = edges[idx]
        if edge[1] == nodeA and edge[2] == nodeB:
            return (edge, idx)
        elif edge[1] == nodeB and edge[2] == nodeA:
            return (edge, idx)
    return None

# TODO maybe don't need SequenceMap
# mutable struct SequenceMap
#     original::String
#     original_index::Int
#     aligned::String
#     aligned_index::Int
# end

# mutable struct Particle
#     sequences::Array{String,1}              # initial order/position
#     aligned_sequences::Array{String,1}      # hold the aligned sequences to avoid unnecessary recomputes
#     #sequence_map::Array{SequenceMap,1}
#     position::Array{String,1} # ordering    # current order/position
#     best_sequences::Array{String,1}         # best order/position
#     best_score::Int                         # best score
#     pidx::Int                               # particle index (for reverse mapping into the particle array)
#     edges::Array{Tuple{String,String,Int},1}# edges for the current position
# end
class Particle:
    def __init__(self):
        self.sequences = []            # initial order/position
        self.aligned_sequences = []    # hold the aligned sequences to avoid unnecessary recomputes
        self.position = [] # ordering  # current order/position
        self.best_sequences = []       # best order/position
        self.best_score = 0            # best score
        self.pidx = 0                  # particle index (for reverse mapping into the particle array)
        self.edges = []                # edges for the current position

def print_particle(p:Particle):
    print("Particle: ", p.pidx)
    print("  position:")
    for seq in p.position:
        print("    ", seq)

    print("  best_sequences:")
    for seq in p.best_sequences:
        print("    ", seq)

    print("  best_score: ", p.best_score)
    print("  edges:")
    for edge in p.edges:
        print("    ", str(edge))

def argmax(l):
    if len(l) == 0:
        raise RuntimeError("called argmax on empty list")
    midx = 0
    for i in range(1, len(l)):
        if l[i] > l[midx]:
            midx = i
    return midx


def PSO_MSA(sequences:list, iterations:int):
    # N = 40 # length of the generated sequences
    # t = 10 # number of sequences
    # iterations = 2000
    t = len(sequences)
    solution_space = math.factorial(t) #/ (t*(t-1))

    num_particles = 0
    if t > 10:
        num_particles = 100
    else:
        num_particles = t*(t-1) # num edges

    print("num_particles: %d" % num_particles)
    search_space = num_particles * iterations
    search_solution_ratio = search_space / solution_space

    print("Solution Space: %.02f" % solution_space)
    print("Search Space: %.02f" % search_space)
    print("Search/Solution space: %.02f" % search_solution_ratio)

    nodes = []
    for i in range(0, len(sequences)):
        nodes.append(sequences[i])

    permutations = random_permutations(nodes, num_particles)

    # Initialize particles list from the generated permutations
    print("Initializing particles")
    particles = [None] * len(permutations)
    for pidx in range(0, len(permutations)):
        sequences = copy.copy(permutations[pidx])
        aligned_sequences = copy.copy(sequences)

        p_edges = []
        for sidx in range(1, len(sequences)):
            seqA = sequences[sidx-1]
            seqB = sequences[sidx]
            scoreAB = score_sequences([seqA, seqB])
            p_edges.append((seqA, seqB, scoreAB))

        aligned_sequences = progressive_alignment_inorder(sequences, p_edges)
        initial_best_score = score_sequences(aligned_sequences)
        p = Particle()

        p.sequences = permutations[pidx]            # original sequences/position
        p.aligned_sequences = aligned_sequences     # aligned_sequences
        p.position = copy.copy(permutations[pidx])       # current_position
        p.best_sequences = copy.copy(permutations[pidx]) # best_sequences
        p.best_score = initial_best_score           # initial best score
        p.pidx = pidx                               # index
        p.edges = p_edges                           # edges for the current position
        particles[pidx] = p

    if debug:
        print("INITIAL PARTICLES")
        for p in particles:
            print_particle(p)

    #alpha = Random.rand()
    #beta = Random.rand()
    alpha = 0.2
    beta = 0.2

    # keep track of the local best score+ordering for each particle
    local_best_scores = [p.best_score for p in particles]
    # keep track of the global best particle
    #best_score = max(local_best_scores)
    global_best_particle = particles[argmax(local_best_scores)]
    print("Initial Best")
    print("Global best (score = %d, pidx = %d):\n" % (global_best_particle.best_score, global_best_particle.pidx))
    for seq in global_best_particle.sequences:
        print(seq)

    print("Initial global best score: %d" % global_best_particle.best_score)

    iteration_metrics = []

    random_subgroup_move_max = 10
    random_subgroup_move_count = 0
    random_subgroup_prob = 0.1

    # Percentage of particles which if they don't move, move them randomly
    reshuffle_threshold = len(particles) * 0.75

    for iteration in range(iterations):
        # Revert all particle aligned sequences to their original unaligned values
        # for p in particles
        #     p.aligned_sequences = copy(p.position)
        # end
        unmoved_particles = []

        print("Beginning iteration %d" % iteration)
        start = datetime.datetime.now()

        if random_subgroup_move_count >= random_subgroup_move_max:
            do_random_subgroup_move = True
            random_subgroup_move_count = 0
        else:
            random_subgroup_move_count += 1
            do_random_subgroup_move = False

        for pidx in range(0, len(particles)):
            particle = particles[pidx]
            #println("Looking at particle:")
            #print_particle(particles[pidx])

            new_velocity = None

            if do_random_subgroup_move:
                if random.random() < random_subgroup_prob:
                    random_position = random_permutations(particle.position, 1)[0]
                    dprintf("Randomly moving particle from %s to %s\n", str(particle.position), str(random_position))
                    new_velocity = pso_particle_velocity(
                        random_position,
                        global_best_particle.best_sequences,
                        particle.best_sequences,
                        alpha, beta)

            if new_velocity == None:
                new_velocity = pso_particle_velocity(
                    particle.position,
                    global_best_particle.best_sequences,
                    particle.best_sequences,
                    alpha, beta)

            if len(new_velocity) == 0:
                unmoved_particles.append(pidx)

            if len(new_velocity) > 0:
                # Apply the velocity to the position of the particle
                new_position = apply_velocity_to_particle(particle.position, new_velocity)
                dprintf("Applied velocity %s to particle %d, moved from\n%s, to %s\n",
                    str(new_velocity),
                    particle.pidx,
                    particle.position,
                    new_position)
                particle.position = new_position

                # Calculate the new edges for the current position (sequence ordering)
                p_edges = [None] * (len(particle.position)-1)
                for sidx in range(1, len(particle.position)):
                    seqA = particle.position[sidx-1]
                    seqB = particle.position[sidx]
                    scoreAB = score_sequences([seqA, seqB])
                    p_edges[sidx-1] = (seqA, seqB, scoreAB)

                particle.edges = p_edges

                # Perform a progressive alignment of the new position with the new edges
                aligned_sequences = progressive_alignment_inorder(particle.position, particle.edges)
                score = score_sequences(aligned_sequences)

                particle.aligned_sequences = aligned_sequences

                # if the new position score is better, update the local best
                if score > particle.best_score:
                    dprintf("Updating particle best score to %d\n", score)
                    particle.best_score = score
                    particle.best_sequences = copy.copy(particle.position)

                # if particle.best_score > global_best_particle.best_score
                #     @printf("Updating global best score to %d\n", particle.best_score)
                #     global_best_particle = deepcopy(particle)
                # end

            # Update the global best at the end of the round
            particle_best_scores = [p.best_score for p in particles]
            best_particle = particles[argmax(particle_best_scores)]
            if best_particle.best_score > global_best_particle.best_score:
                print("Updating global best score to %d\n" % best_particle.best_score)
                global_best_particle = copy.deepcopy(best_particle)

        if len(unmoved_particles) >= reshuffle_threshold:
            dprintf("Reshuffling particles: %s\n", str(unmoved_particles))
            for reshuffle_pidx in unmoved_particles:
                p = particles[reshuffle_pidx]
                p.position = random_permutations(p.position, 1)[0]

        end_time = datetime.datetime.now()
        duration = end_time - start
        iteration_metrics.append(duration)
        # @printf("\nIteration %d took %d milliseconds\n", iteration, duration.value)
        # @printf("Global best (score = %d, pidx = %d):\n",
        #     global_best_particle.best_score,
        #     global_best_particle.pidx)
        # for seq in global_best_particle.aligned_sequences
        #     println(seq)
        # end
        # println()

    print("FINISHED PARTICLE SWARM")

    # print the final results
    print("\nGlobal best (score = %d, pidx = %d):\n" % (
        global_best_particle.best_score,
        global_best_particle.pidx))
    for seq in global_best_particle.aligned_sequences:
        print(seq)
    print()

    return global_best_particle.best_score

def main():
    t = 5 # number of sequences
    N = 38
    sequences = generate_sequences(t, N)
    print(sequences)
    perf = PSO_MSA(sequences, 200)
    #perf = @timed PSO_MSA(sequences, 200)
    #println(perf)

main()