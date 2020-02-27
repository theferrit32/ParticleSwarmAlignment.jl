
using Printf

include("./ParticleSwarmAlignment.jl")

function main()
    #GC.enable(false)
    if length(ARGS) > 0
        if length(ARGS) != 3
            println("Usage:")
            println("  julia main.jl number_of_sequences length_of_sequences swarm_iterations")
            throw(ErrorException("Command line args were provided but there were not 3 of them"))
        end
        t = parse(Int64, ARGS[1])
        N = parse(Int64, ARGS[2])
        iterations = parse(Int64, ARGS[3])
        @printf("Using provided args t=%d, N=%d, iterations=%d\n", t, N, iterations)
    else
        # number of sequences
        t = 5
        # length of sequences
        N = 38
        # iterations
        iterations = 200
        @printf("Using default args t=%d, N=%d, iterations=%d\n", t, N, iterations)
    end

    sequences = ParticleSwarmAlignment.generate_sequences(t, N)
    println(sequences)
    val, t, bytes, gctime, memallocs = @timed ParticleSwarmAlignment.PSO_MSA(sequences, iterations)

    @printf("Total bytes: %d\n", bytes)
end

main()
