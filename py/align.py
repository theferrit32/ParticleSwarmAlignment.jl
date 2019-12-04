
def argmax(l):
    if len(l) == 0:
        raise RuntimeError("called argmax on empty list")
    midx = 0
    for i in range(1, len(l)):
        if l[i] > l[midx]:
            midx = i
    return midx

def global_align(v, w, match_penalty=1, mismatch_penalty=-0.5, deletion_penalty=-2):
    #print("global_align, v = %s, w = %s" % (v, w))
    vlen = len(v)
    wlen = len(w)
    s = []
    b = []
    for i in range(wlen+1):
        s_row = []
        b_row = []
        for j in range(vlen+1):
            s_row.append(-1)
            b_row.append(-1)
        s.append(s_row)
        b.append(b_row)
    assert(len(s) == wlen+1)
    assert(len(s[0]) == vlen+1)

    # matrix is v along rows, w along columns
    #                 v (n1+1)
    #          [                    ]
    #          [                    ]
    #  w(n2+1) [                    ]
    #          [                    ]

    # b values
    # diag = 0
    # up   = 1
    # left = 2

    # fill first row
    for j in range(0, vlen+1):
        s[0][j] = (j) * deletion_penalty
        b[0][j] = 2

    # fill first column
    for i in range(0, wlen+1):
        s[i][0] = (i) * deletion_penalty
        b[i][0] = 1

    # for b_row in b:
    #     print(b_row)

    for i in range(1, wlen+1):
        for j in range(1, vlen+1):
            #print("i =", i, "j =", j)
            if v[j-1] == w[i-1]:
                ms = s[i-1][j-1] + match_penalty
            else:
                ms = s[i-1][j-1] + mismatch_penalty

            test = [ms, s[i-1][j] + deletion_penalty, s[i][j-1] + deletion_penalty]
            p = argmax(test)
            s[i][j] = test[p]
            b[i][j] = p

    i = wlen # row idx
    j = vlen # col idx
    sv = []
    sw = []
    # for b_row in b:
    #     print(b_row)
    #print("Reconstructing alignment")
    while(i > 0 or j > 0):
        #print("i =", i, "j =", j)
        p = b[i][j]
        if (p == 0): # diag
            i = i-1
            j = j-1
            sv.append(v[j])
            sw.append(w[i])
        elif p == 1: # up
            i = i-1
            sv.append("-")
            sw.append(w[i])
            #sv.append(v[i])
            #sw.append("-")
        elif p == 2: # left
            j = j-1
            sv.append(v[j])
            sw.append("-")
            #sv.append("-")
            #sw.append(w[j])
        else:
            break

    sv.reverse()
    sw.reverse()
    v_out = "".join(sv)
    w_out = "".join(sw)
    #print("sv =", sv)
    #print("sw =", sw)
    #print("v_out =", v_out)
    #print("w_out =", w_out)
    return (s[wlen][vlen], v_out, w_out)
