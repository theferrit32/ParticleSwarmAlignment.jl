#!/bin/bash
# This script uses /usr/bin/time to obtain the max-resident memory for each process
set -e -x

# NOTE /usr/bin/time memory values are kilobytes
# Vary the length of the input sequences between 20 and 500
# syntax for range: seq first increment last

# Number of sequences to align
num_sequences=5
# Number of particle swarm iterations
iterations=3
# Number of tests to perform for each parameter value
test_iters=10

run python test
rm -f python-mem.txt python-time.txt python-output.txt
for length in `seq 20 20 100`; do
    for i in `seq 1 $test_iters`; do
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
    done
done


# run julia test
rm -f julia-mem.txt julia-time.txt julia-output.txt
for length in `seq 20 20 100`; do
    for i in `seq 1 $test_iters`; do
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
    done
done