using Dates
start_time = Dates.now()
println("Starting imports: ", start_time)
# using StatsBase
# using Distributions
# using Gadfly
# using TimerOutputs
using Random
#import Cairo, Fontconfig
using JSON
end_time = Dates.now()
println("Finished imports: ", end_time - start_time)
include("./main.jl")

# Function to plot graph
# function plotgraph(xaxis::Array{Int64,1},yaxis::Array{Float64,1},x_label::String,y_label::String,heading::String)
#     e = plot(x=xaxis, y=yaxis, Guide.xlabel(x_label), Guide.ylabel(y_label), Guide.title(heading), Geom.line,
#         Theme(minor_label_font_size=14pt, major_label_font_size=14pt));
#     return(e)
# end

function write_results(filename, x_vals, scores, mem, time)
    f = open(filename, write=true)
    d = Dict(
        "x_vals" => x_vals,
        "scores" => scores,
        "mems" => mem,
        "times" => time
    )
    JSON.write(f, json(d))
end

#function for testing with increasing number of input sequences.
function msatest1()
    println("TEST msatest1")
    score_mat = Float64[]
    num_of_inp_seq = Int64[]
    bytes = Float64[]
    time = Float64[]

    #Vary the number of input sequences and record reading.
    for t in 10:20:500
        println("t=", t)
        sequences = generate_sequences(t, 50)
        MSA = @timed PSO_MSA(sequences,20)
        push!(score_mat,MSA[1])
        push!(num_of_inp_seq,t)
        push!(bytes, MSA[3])
        push!(time, MSA[2] - MSA[4])
    end

    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(num_of_inp_seq)
    popfirst!(bytes)
    popfirst!(time)

    println(num_of_inp_seq)
    println(score_mat)
    println(bytes)

    write_results("msatest1.json", num_of_inp_seq, score_mat, bytes, time)

    #z1 = plotgraph(num_of_inp_seq, score_mat, "Number of input sequences", "Score", "Plot of Score v/s Number of input sequences")
    #z2 = plotgraph(num_of_inp_seq, bytes, "Number of input sequences", "Memory in bytes", "Plot of Memory v/s Number of input sequences")
    #z3 = plotgraph(num_of_inp_seq, time, "Number of input sequences", "Time in seconds", "Plot of Time v/s Number of input sequences")
    #display(z1)
    #display(z2)
    #display(z3)

    #draw(PNG("msatest1-z1.png", 8inch, 8inch), z1)
end

#function for testing with increasing length of input sequences.
function msatest2()
    println("TEST msatest2")
    score_mat = Float64[]
    inp_seq_len = Int64[]
    bytes = Float64[]
    time = Float64[]

    #Vary the length of input sequences and record reading.
    for N in 20:20:800
        println("N=", N)
        sequences = generate_sequences(10, N)
        MSA = @timed PSO_MSA(sequences, 20)
        push!(score_mat, MSA[1])
        push!(inp_seq_len, N)
        push!(bytes, MSA[3])
        push!(time, MSA[2]-MSA[4])
    end

    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(inp_seq_len)
    popfirst!(bytes)
    popfirst!(time)

    write_results("msatest2.json", inp_seq_len, score_mat, bytes, time)

    # z1 = plotgraph(inp_seq_len, score_mat, "Input sequence length", "Score", "Plot of Score v/s Length of input sequence")
    # z2 = plotgraph(inp_seq_len, bytes, "Input sequence length", "Memory in bytes", "Plot of Memory v/s Length of input sequence")
    # z3 = plotgraph(inp_seq_len, time, "Input sequence length", "Time in seconds", "Plot of Time v/s Length of input sequence")
    # display(z1)
    # display(z2)
    # display(z3)
end

#function for testing with increasing number of iterations.
function msatest3()
    println("TEST msatest3")
    score_mat = Float64[]
    itr = Int64[]
    bytes = Float64[]
    time = Float64[]
    sequences = generate_sequences(10, 50)

    #Vary the number of iterations and record reading.
    for i in 0:5:120
        println("i=", i)
        Random.seed!(1)
        MSA = @timed PSO_MSA(sequences,i)
        push!(score_mat,MSA[1])
        push!(itr,i)
        push!(bytes, MSA[3])
        push!(time, MSA[2]-MSA[4])
    end

    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(itr)
    popfirst!(bytes)
    popfirst!(time)

    write_results("msatest3.json", itr, score_mat, bytes, time)

    # z1 = plotgraph(itr, score_mat, "Number of iterations", "Score", "Plot of Score v/s Number of iterations")
    # z2 = plotgraph(itr, bytes, "Number of iterations", "Memory in bytes", "Plot of Memory v/s Number of iterations")
    # z3 = plotgraph(itr, time, "Number of iterations", "Time in seconds", "Plot of Time v/s Number of iterations")
    # display(z1)
    # display(z2)
    # display(z3)
end

#msatest1()
#msatest2()
msatest3()
