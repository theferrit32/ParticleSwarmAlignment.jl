#!/bin/bash
# This script uses /usr/bin/time to obtain the max-resident memory for each process
set -e -x

# NOTE memory values are kilobytes
# We vary the length of the input sequences between 20 and 500
# syntax for range: seq first increment last
# Number of sequences to align
num_sequences=5
# Number of particle swarm iterations
iterations=3
# Number of tests to perform for each parameter value
test_iters=5

# run python test
rm -f python-mem.txt python-time.txt python-output.txt
for length in `seq 20 20 400`; do
    time_counter=0
    internal_bytes_counter=0
    max_res_counter=0

    for i in `seq 1 $test_iters`; do
        echo "LENGTH $length"
        start=$(date +%s.%N)
        output=$(/usr/bin/time -f "%M" python py/main.py $num_sequences $length $iterations 2>&1)
        end=$(date +%s.%N)
        duration=$(echo "$end - $start" | bc -l)
        printf "LENGTH $length\n" >> python-output.txt
        printf "$output\n" >> python-output.txt
        internal_bytes=$(echo "$output" | tail -n 2 | head -n 1 | awk '{print $3}') # second to last line, 3rd term
        max_res=$(echo "$output" | tail -n 1)
        printf "$length $internal_bytes $max_res\n" >> python-mem.txt
        printf "$length $duration\n" >> python-time.txt

        # time_counter=$(echo "$time_counter + $duration" | bc -l)
        # internal_bytes_counter=$(echo "$internal_bytes_counter + $internal_bytes" | bc -l)
        # max_res_counter=$(echo "$max_res_counter + $max_res" | bc -l)
    done
    # avg_time=$(echo "$time_counter / $iters" | bc -l)
    # avg_internal_bytes=$(echo "$internal_bytes_counter / $iters" | bc -l)
    # avg_max_res=$(echo "$max_res_counter / $iters" | bc -l)
    # printf "$length $avg_internal_bytes $avg_max_res\n" >> python-mem.txt
    # printf "$length $avg_time\n" >> python-time.txt
done


# run julia test
rm -f julia-mem.txt julia-time.txt julia-output.txt
for length in `seq 20 20 400`; do
    time_counter=0
    internal_bytes_counter=0
    max_res_counter=0

    for i in `seq 1 $test_iters`; do
        echo "LENGTH $length"
        start=$(date +%s.%N)
        output=$(/usr/bin/time -f "%M" julia main.jl $num_sequences $length $iterations 2>&1)
        end=$(date +%s.%N)
        duration=$(echo "$end - $start" | bc -l)
        printf "LENGTH $length\n" >> julia-output.txt
        printf "$output\n" >> julia-output.txt
        internal_bytes=$(echo "$output" | tail -n 2 | head -n 1 | awk '{print $3}') # second to last line, 3rd term
        max_res=$(echo "$output" | tail -n 1)
        printf "$length $internal_bytes $max_res\n" >> julia-mem.txt
        printf "$length $duration\n" >> julia-time.txt

        # time_counter=$(echo "$time_counter + $duration" | bc -l)
        # internal_bytes_counter=$(echo "$internal_bytes_counter + $internal_bytes" | bc -l)
        # max_res_counter=$(echo "$max_res_counter + $max_res" | bc -l)
    done
    # avg_time=$(echo "$time_counter / $iters" | bc -l)
    # avg_internal_bytes=$(echo "$internal_bytes_counter / $iters" | bc -l)
    # avg_max_res=$(echo "$max_res_counter / $iters" | bc -l)
    # printf "$length $avg_internal_bytes $avg_max_res\n" >> julia-mem.txt
    # printf "$length $avg_time\n" >> julia-time.txt
done