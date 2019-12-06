import json
import matplotlib.pyplot as plt
font = {'family' : 'sans',
        'size'   : 18}
plt.rc('font', **font)

#plt.rcParams['figure.figsize'] = 8, 6
#plt.margins(x=0.005, y=0.005, tight=True)

files = [
    'msatest1.json',
    'msatest2.json',
    'msatest3.json']

def test1():
    with open('msatest1.json') as f:
        msatest1_obj = json.load(f)

    msatest1_obj['mems'] = [i/(1024*1024*1024) for i in msatest1_obj['mems']]

    fig, ax = plt.subplots(1, 3, squeeze=True, figsize=figsize)

    ax[0].plot(msatest1_obj['x_vals'], msatest1_obj['scores'], marker='o')
    ax[1].plot(msatest1_obj['x_vals'], msatest1_obj['mems'], marker='o')
    ax[2].plot(msatest1_obj['x_vals'], msatest1_obj['times'], marker='o')
    ax[0].set_ylabel('Score', labelpad=ylabelpad)
    ax[1].set_ylabel('Cumulative memory (GiB)', labelpad=ylabelpad)
    ax[2].set_ylabel('Time (seconds)', labelpad=ylabelpad)
    ax[0].set_xlabel('Number of Input Sequences')
    ax[1].set_xlabel('Number of Input Sequences')
    ax[2].set_xlabel('Number of Input Sequences')
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()

    plt.tight_layout()
    plt.show()
def test2():
    with open('msatest2.json') as f:
        msatest2_obj = json.load(f)

    msatest2_obj['mems'] = [i/(1024*1024*1024) for i in msatest2_obj['mems']]

    fig, ax = plt.subplots(1, 3, squeeze=True, figsize=figsize)

    ax[0].plot(msatest2_obj['x_vals'], msatest2_obj['scores'], marker='o')
    ax[1].plot(msatest2_obj['x_vals'], msatest2_obj['mems'], marker='o')
    ax[2].plot(msatest2_obj['x_vals'], msatest2_obj['times'], marker='o')
    ax[0].set_ylabel('Score', labelpad=ylabelpad)
    ax[1].set_ylabel('Cumulative memory (GiB)', labelpad=ylabelpad)
    ax[2].set_ylabel('Time (seconds)', labelpad=ylabelpad)
    ax[0].set_xlabel('Length of Input Sequences')
    ax[1].set_xlabel('Length of Input Sequences')
    ax[2].set_xlabel('Length of Input Sequences')
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()

    plt.tight_layout()
    plt.show()
def test3():
    with open('msatest3.json') as f:
        msatest3_obj = json.load(f)

    msatest3_obj['mems'] = [i/(1024*1024*1024) for i in msatest3_obj['mems']]

    fig, ax = plt.subplots(1, 3, squeeze=True, figsize=figsize)

    ax[0].plot(msatest3_obj['x_vals'], msatest3_obj['scores'], marker='o')
    ax[2].plot(msatest3_obj['x_vals'], msatest3_obj['mems'], marker='o')
    ax[1].plot(msatest3_obj['x_vals'], msatest3_obj['times'], marker='o')

    ax[0].set_ylabel('Score', labelpad=ylabelpad)
    ax[2].set_ylabel('Cumulative memory (GiB)', labelpad=ylabelpad)
    ax[1].set_ylabel('Time (seconds)', labelpad=ylabelpad)

    ax[0].set_xlabel('Number of PSO Iterations')
    ax[2].set_xlabel('Number of PSO Iterations')
    ax[1].set_xlabel('Number of PSO Iterations')

    ax[0].grid()
    ax[2].grid()
    ax[1].grid()

    plt.tight_layout()
    plt.show()

ylabelpad = 10
figsize=(18, 4.75)
test1()
test2()
test3()