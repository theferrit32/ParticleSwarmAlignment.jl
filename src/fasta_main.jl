include("./ParticleSwarmAlignment.jl")
using .ParticleSwarmAlignment
using Combinatorics
using Printf

function load_fasta_file(fname)
    sequences = String[]
    line_counter = 0
    open(fname) do f
        seq = ""
        for line in eachline(f)
            if startswith(line, ">")
                if length(seq) > 0
                    push!(sequences, seq)
                end
                seq = ""
            else
                seq = string(seq, line)
            end
        end
        if length(seq) > 0
            push!(sequences, seq)
        end
    end
    return sequences
end

function bigmain()
    # Load single-sequence FASTA files for different species
    files = [
        "california-sea-lion.fasta",
        "cattle.fasta",
        "chicken.fasta",
        "chimpanzee.fasta",
        #"fission-yeast.fasta",
        "human.fasta",
        "long-finned-pilot-whale.fasta",
        "norway-rat.fasta",
        "pacific-white-sided-dolphin.fasta",
        "rhesus-monkey.fasta",
        "tropical-clawed-frog.fasta",
        "water-buffalo.fasta",
        "white-tufted-ear-marmoset.fasta"
    ]
    for i in 1:length(files)
        files[i] = string("../fasta/HBA1/", files[i])
    end
    sequences = String[]

    for fname in files
        append!(sequences, load_fasta_file(fname))
    end

    iterations = 50
    PSO_MSA(sequences, iterations)
end

function smallmain()
    filename = "../fasta/ENA-ebi.fasta"
    sequences = load_fasta_file(filename)
    iterations = 50
    score, aligned_sequences = PSO_MSA(sequences, iterations)

    println("Performing exhaustive search")
    perms = permutations(sequences)
    total_perms = length(perms)
    best_score = 0
    best_perm = String[]
    best_aligned = String[]
    i = 1
    for perm in perms
        @printf("Exhaustive Search step %d/%d\n", i, total_perms)
        p_edges = Array{Tuple{String,String,Int},1}(undef,length(perm)-1)
        for sidx in 2:length(perm)
            seqA_orig = perm[sidx-1]
            seqB_orig = perm[sidx]
            align_score, seqA_aligned, seqB_aligned = global_align(seqA_orig, seqB_orig)
            scoreAB = score_sequences([seqA_aligned, seqB_aligned])
            p_edges[sidx-1] = (seqA_orig, seqB_orig, scoreAB)
        end
        aligned_sequences = progressive_alignment_inorder(perm, p_edges)
        score = score_sequences(aligned_sequences)
        if score > best_score
            best_score = score
            best_perm = copy(perm)
            best_aligned = copy(aligned_sequences)
        end
        i += 1
    end
    @printf("BEST SEQUENCES (score=%d)\n", best_score)
    for seq in best_aligned
        println(seq)
    end
end

function scoremain()
    sequences = [
        "-------AT-T----TAA-CT--CAGA---ACTGTC-GG---T--A-GACC----G-TG-A--AG--G-CGT--TCCC",
        "--------T-T----T---CC--CAGAG--AG-GTA-GGCTTT--ATGAGC----G--G-ATTAG--G--GTGAT---",
        "---T-ACAC-CTTAGTG--CT--CTCA---AT-GT--GG--TTG-A-GATCT---G-TT-A--A-----C-T------",
        "-AGTGA-AC-C---GTG--CG--C--A---AT--T------TTGTA-CATCT---G-TA-AG-AG--G-C-TG-T---",
        "-GGAGG-AG-C----AGA-CT--C------AC--TA-GG---T--G-TAT-----G----AG-AG--G-C-TG-T---",
        "---C-G-GTGT----TG-TCCGTC-GGC--AC-CTG-AG---TG-A---CTT---GATGCAT-------C-T--T---",
        "TGAC-G-GC-T----TGATCC--C-G----A---TG-GG--CTG-A--ACC----GATACAT-AG--G-C-C------",
        "-----A-AC-CC---AG---T--C-AA---AC-ATACGGA-ACC-A-CGCCTACAG-C--A--AG--G-C----C---",
        "-------AC------AGG-GT--CTAA---G--GTAC----TCC-T-TGTCCAC-G-C--A--TGACGTC-TA-C---",
        "-----ATGC-C-----GC-TT--C-GGTTTA--GTAC----TTG-A--GTCC-CAG-CG-A--AG-CG-C--G-----"
    ]
    score = score_sequences(sequences)
    println(score)
end

#scoremain()
#bigmain()
smallmain()