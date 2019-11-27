using StatsBase
using Distributions
using Gadfly
using TimerOutputs
include("./main.jl")

# Function to plot graph
function plotgraph(xaxis::Array{Int64,1},yaxis::Array{Float64,1},x_label::String,y_label::String,heading::String)
    e = plot(x=xaxis, y=yaxis, Guide.xlabel(x_label), Guide.ylabel(y_label), Guide.title(heading), Geom.line)
    #plot(e)
    return(e)
end

#function for testing with increasing number of input sequences.
function msatest1()
    score_mat = Float64[]
    num_of_inp_seq = Int64[]
    bytes = Float64[]
    time = Float64[]
    
    #Vary the length of input sequences and record reading.
    for t in 10:5:50
        sequences = generate_sequences(t, 50)
        MSA = @timed PSO_MSA(50,t,1,sequences)
        push!(score_mat,MSA[1]/t)
        push!(num_of_inp_seq,t)
        push!(bytes, MSA[3]/t)
        push!(time, MSA[2] - MSA[4])
    end
    
    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(num_of_inp_seq)
    popfirst!(bytes)
    popfirst!(time)
    
    z1 = plotgraph(num_of_inp_seq, score_mat, "Number of input sequences", "Normalized Score", "Plot of Normalized Score v/s Number of input sequences")
    z2 = plotgraph(num_of_inp_seq, bytes, "Number of input sequences", "Normalized space", "Plot of Normalized Space v/s Number of input sequences")
    z3 = plotgraph(num_of_inp_seq, time, "Number of input sequences", "Time", "Plot of Time v/s Number of input sequences")

    display(z1)
    display(z2)
    display(z3)
    
end
#function for testing with increasing length of input sequences.
function msatest2()
    score_mat = Float64[]
    inp_seq_len = Int64[]
    bytes = Float64[]
    time = Float64[]
    
    #Vary the length of input sequences and record reading.
    for N in 50:10:100
        sequences = generate_sequences(10, N)
        MSA = @timed PSO_MSA(N,10,1,sequences)
        push!(score_mat,MSA[1]/N)
        push!(inp_seq_len,N)
        push!(bytes, MSA[3]/N)
        push!(time, MSA[2]-MSA[4])
    end
    
    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(inp_seq_len)
    popfirst!(bytes)
    popfirst!(time)
    
    z1 = plotgraph(inp_seq_len, score_mat, "Input sequence length", "Normalized Score", "Plot of Normalized Score v/s Length of input sequence")
    z2 = plotgraph(inp_seq_len, bytes, "Input sequence length", "Normalized space", "Plot of Normalized Space v/s Length of input sequence")
    z3 = plotgraph(inp_seq_len, time, "Input sequence length", "Time", "Plot of Time v/s Length of input sequence")

    display(z1)
    display(z2)
    display(z3)

end
#function for testing with increasing number of iterations.
function msatest3()
    score_mat = Float64[]
    itr = Int64[]
    bytes = Float64[]
    time = Float64[]
    sequences = generate_sequences(5, 50)
    
    #Vary the length of input sequences and record reading.
    for i in 1:2:13
        MSA = @timed PSO_MSA(50,5,i,sequences)
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
    
    z1 = plotgraph(itr, score_mat, "Number of iterations", "Score", "Plot of Score v/s Number of iterations")
    z2 = plotgraph(itr, bytes, "Number of iterations", "Space", "Plot of Space v/s Number of iterations")
    z3 = plotgraph(itr, time, "Number of iterations", "Time", "Plot of Time v/s Number of iterations")

    display(z1)
    display(z2)
    display(z3)

end

msatest1()
msatest2()
msatest3()