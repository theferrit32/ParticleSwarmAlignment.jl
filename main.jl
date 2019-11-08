using Random
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

sequences = generate_sequences(5, 10)
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
            #function global_align(v,w,match_penalty=1,mismatch_penalty=-1,deletion_penalty=0)
            distanceAB, seqA_aligned, seqB_aligned = global_align(seqA, seqB, 1, -1, 0)
            push!(edges, (seqA, seqB, distanceAB))
        end
    end
end

print("Nodes: ")
println(nodes)
print("Edges: ")
println(edges)