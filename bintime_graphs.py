import json
import statistics
import math
import matplotlib.pyplot as plt
font = {'family' : 'sans',
        'size'   : 14}
plt.rc('font', **font)
#plt.rcParams['legend.loc'] = 'upper left'

def group_and_stdev(xinputs, yinputs):
    assert(len(xinputs) == len(yinputs))

    groupings = {}
    for x,y in zip(xinputs, yinputs):
        if x not in groupings:
            groupings[x] = []
        groupings[x].append(y)

    x_vals = []
    y_vals = []
    intervals = []

    # find the standard deviation for each group
    for x,ys in groupings.items():
        n = len(ys)
        mean = sum(ys) / len(ys)
        sdev = statistics.stdev(ys, mean)
        z = 1.960 # 95% interval
        interval = z * sdev / (math.sqrt(n))

        # store the results
        x_vals.append(x)
        y_vals.append(mean)
        intervals.append(interval)

    return x_vals, y_vals, intervals

def load_mem_file(fname):
    x_vals = []
    y_totals = []
    y_maxes = []

    with open(fname) as fin:
        for line in fin:
            line_terms = line.strip().split(" ")
            sequence_lengths = int(line_terms[0])
            total_bytes = int(line_terms[1])
            max_resident = int(line_terms[2])
            x_vals.append(sequence_lengths)
            y_totals.append(total_bytes)
            y_maxes.append(max_resident)

    x_vals_out, y_totals, totals_intervals = group_and_stdev(x_vals, y_totals)
    x_vals_out, y_maxes, maxes_intervals = group_and_stdev(x_vals, y_maxes)

    return x_vals_out, y_totals, totals_intervals, y_maxes, maxes_intervals

def load_time_file(fname):
    x_vals = []
    y_times = []
    with open(fname) as fin:
        for line in fin:
            line_terms = line.strip().split(" ")
            sequence_lengths = int(line_terms[0])
            t = float(line_terms[1])
            x_vals.append(sequence_lengths)
            y_times.append(t)

    # group them by the x_val
    groupings = {}
    for x,y in zip(x_vals, y_times):
        if x not in groupings:
            groupings[x] = []
        groupings[x].append(y)

    x_vals = []
    y_vals = []
    intervals = []

    # find the standard deviation for each group
    for x,ys in groupings.items():
        n = len(ys)
        mean = sum(ys) / len(ys)
        sdev = statistics.stdev(ys, mean)
        z = 1.960 # 95% interval
        interval = z * sdev / (math.sqrt(n))

        # store the results
        x_vals.append(x)
        y_vals.append(mean)
        intervals.append(interval)

    return x_vals, y_vals, intervals

def mem_plots_total():
    fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(12, 6))

    julia_x_vals, julia_y_totals, julia_total_ivls, julia_y_maxes, julia_max_ivls = load_mem_file("julia-mem.txt")
    print("Dividing julia internal memory usage by 1024 to get kilobytes")
    julia_y_totals = [i / 1024 for i in julia_y_totals]
    #ax.plot(julia_x_vals, julia_y_totals, label="Total Memory Allocated (Julia)")
    ax.errorbar(julia_x_vals, julia_y_totals,
        yerr=julia_total_ivls,
        label="Total Memory Allocated (Julia)",
        color=color1,
        ecolor=ecolor1)

    python_x_vals, python_y_totals, python_total_ivls, python_y_maxes, python_max_ivls = load_mem_file("python-mem.txt")
    #ax.plot(python_x_vals, python_y_totals, label="Total Memory Allocated (Python)")
    ax.errorbar(python_x_vals, python_y_totals,
        yerr=python_total_ivls,
        label="Total Memory Allocated (Python)",
        color=color2,
        ecolor='red')

    ax.set_xlabel("Input Sequence Lengths")
    ax.set_ylabel("Memory (KiB)")
    #ax.set_ylim((0,1500000))
    ax.grid()
    plt.legend()
    plt.title("Total Memory Allocated by PSO-MSA Algorithm")
    plt.tight_layout()
    plt.show()

def mem_plots_max():
    fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(12, 6))

    julia_x_vals, julia_y_totals, julia_total_ivls, julia_y_maxes, julia_max_ivls = load_mem_file("julia-mem.txt")
    print("Dividing julia internal memory usage by 1024 to get kilobytes")
    julia_y_totals = [i / 1024 for i in julia_y_totals]
    #ax.plot(julia_x_vals, julia_y_maxes, label="Maximum Process Memory Footprint (Julia)")
    ax.errorbar(julia_x_vals, julia_y_maxes,
        yerr=julia_max_ivls,
        label="Maximum Process Memory Footprint (Julia)",
        color=color1,
        ecolor='red')

    python_x_vals, python_y_totals, python_total_ivls, python_y_maxes, python_max_ivls = load_mem_file("python-mem.txt")
    #ax.plot(python_x_vals, python_y_maxes, label="Maximum Process Memory Footprint (Python)")
    ax.errorbar(python_x_vals, python_y_maxes,
        yerr=python_max_ivls,
        label="Maximum Process Memory Footprint (Python)",
        color=color2,
        ecolor='red')

    ax.set_xlabel("Input Sequence Lengths")
    ax.set_ylabel("Memory (KiB)")
    #ax.set_ylim((0,1500000))
    ax.grid()
    plt.legend()
    plt.title("Maximum Memory Space Allocated to PSO-MSA Process")
    plt.tight_layout()
    plt.show()


def time_plots():
    # must turn off tracemalloc in python to get accurate measure
    fig, ax = plt.subplots(1, 1, squeeze=True, figsize=(12, 6))

    julia_x_vals, julia_y_time, julia_intervals = load_time_file("julia-time.txt")
    # ax.plot(julia_x_vals, julia_y_time, color='yellow', label="Process Execution Time (Julia)")
    ax.errorbar(julia_x_vals, julia_y_time,
        yerr=julia_intervals,
        label="Process Execution Time (Julia)",
        color=color1,
        ecolor=ecolor1)

    python_x_vals, python_y_time, python_intervals = load_time_file("python-time.txt")
    #ax.plot(python_x_vals, python_y_time, label="Process Execution Time (Python)", color='yellow')
    ax.errorbar(python_x_vals, python_y_time,
        yerr=python_intervals,
        label="Process Execution Time (Python)",
        color=color2,
        ecolor=ecolor1)

    ax.set_xlabel("Input Sequence Lengths")
    ax.set_ylabel("Time (Seconds)")
    #ax.set_ylim((0,1500000))
    ax.grid()
    plt.legend()
    plt.title("Execution Time of PSO-MSA Process")
    plt.tight_layout()
    plt.show()


color1 = 'tab:blue'
color2 = 'tab:green'
ecolor1 = 'tab:red'
time_plots()
mem_plots_total()
mem_plots_max()