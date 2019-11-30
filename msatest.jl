using StatsBase
using Distributions
using Gadfly
using TimerOutputs
using Random
include("./main.jl")

# Function to plot graph
function plotgraph(xaxis::Array{Int64,1},yaxis::Array{Float64,1},x_label::String,y_label::String,heading::String)
    e = plot(x=xaxis, y=yaxis, Guide.xlabel(x_label), Guide.ylabel(y_label), Guide.title(heading), 
        Theme(minor_label_font_size=14pt,major_label_font_size=14pt),Geom.line)
    #plot(e)
    return(e)
end

#function for testing with increasing number of input sequences.
function msatest1()
    score_mat = Float64[]
    num_of_inp_seq = Int64[]
    bytes = Float64[]
    time = Float64[]
    
    #Vary the number of input sequences and record reading.
    for t in 10:10:100
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
    
    z1 = plotgraph(num_of_inp_seq, score_mat, "Number of input sequences", "Score", "Plot of Score v/s Number of input sequences")
    z2 = plotgraph(num_of_inp_seq, bytes, "Number of input sequences", "Memory in bytes", "Plot of Memory v/s Number of input sequences")
    z3 = plotgraph(num_of_inp_seq, time, "Number of input sequences", "Time in seconds", "Plot of Time v/s Number of input sequences")

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
    for N in 30:30:300
        print("N=",N)
        sequences = generate_sequences(10, N)
        MSA = @timed PSO_MSA(sequences,20)
        push!(score_mat,MSA[1])
        push!(inp_seq_len,N)
        push!(bytes, MSA[3])
        push!(time, MSA[2]-MSA[4])
    end
    
    #First reading considered as a dry run.
    popfirst!(score_mat)
    popfirst!(inp_seq_len)
    popfirst!(bytes)
    popfirst!(time)
    
    z1 = plotgraph(inp_seq_len, score_mat, "Input sequence length", "Score", "Plot of Score v/s Length of input sequence")
    z2 = plotgraph(inp_seq_len, bytes, "Input sequence length", "Memory in bytes", "Plot of Memory v/s Length of input sequence")
    z3 = plotgraph(inp_seq_len, time, "Input sequence length", "Time in seconds", "Plot of Time v/s Length of input sequence")

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
    sequences = generate_sequences(10, 50)
    
    #Vary the number of iterations and record reading.
    for i in 20:20:120
        print("i =", i)
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
    
    z1 = plotgraph(itr, score_mat, "Number of iterations", "Score", "Plot of Score v/s Number of iterations")
    z2 = plotgraph(itr, bytes, "Number of iterations", "Memory in bytes", "Plot of Memory v/s Number of iterations")
    z3 = plotgraph(itr, time, "Number of iterations", "Time in seconds", "Plot of Time v/s Number of iterations")

    display(z1)
    display(z2)
    display(z3)

end

#msatest1()
#msatest2()
#msatest3()