# global alignment via dynamic programming

# TODO for mismatch_penalty, accept a matrix like BLOSUM62 for scoring the
# mismatch, instead of a static value
function global_align(v, w, match_penalty=1, mismatch_penalty=-0.5, deletion_penalty=-2)
    n1 = length(v)
    n2 = length(w)
    s = zeros(Float64, n1+1, n2+1)
    b = zeros(Float64, n1+1, n2+1)

    for i in 1:(n1+1)
        s[i,1] = (i-1)*deletion_penalty
        b[i,1] = 2
    end
    for j in 1:(n2+1)
        s[1,j] = (j-1)*deletion_penalty
        b[1,j] = 3
    end

    for i in 2:(n1+1)
        for j in 2:(n2+1)
            if v[i-1] == w[j-1]
                ms = s[i-1,j-1] + match_penalty
            else
                # ignore cases where a letter is paired with a gap
                # do not consider this a mismatch
                # if v[i-1] != '-' && w[j-1] != '-'
                #     ms = s[i-1,j-1] + mismatch_penalty
                # else
                #     # if a letter is paired with a gap, add no penalty
                #     ms = s[i-1,j-1] #+ match_penalty # + 0.5 * mismatch_penalty
                # end
                ms = s[i-1,j-1] + mismatch_penalty
            end
            test = [ms, s[i-1,j] + deletion_penalty, s[i,j-1] + deletion_penalty]
            p = argmax(test)
            s[i,j] = test[p]
            b[i,j] = p
        end
    end

    i = n1+1
    j = n2+1
    sv = []
    sw = []
    while(i>1 || j>1)
        p = b[i,j]
        if (p == 1)
            i = i-1
            j = j-1
            push!(sv, v[i])
            push!(sw, w[j])
        elseif p == 2
            i=i-1
            push!(sv, v[i])
            push!(sw, "-")
        elseif p == 3
            j = j-1
            push!(sv, "-")
            push!(sw, w[j])
        else
            break
        end
    end

    return (s[n1+1,n2+1], join(reverse(sv)), join(reverse(sw)))
end


function needleman_wunsch_msa(v, w, match_penalty=1, mismatch_penalty=-1, deletion_penalty=0)
    n1 = length(v)
    n2 = length(w)
    s = zeros(Int64, n1+1, n2+1)
    b = zeros(Int64, n1+1, n2+1)

    for i in 1:(n1+1)
        s[i,1] = (i-1) * deletion_penalty
        b[i,1] = 2
    end
    for j in 1:(n2+1)
        s[1,j] = (j-1) * deletion_penalty
        b[1,j] = 3
    end

    for i in 2:(n1+1)
        for j in 2:(n2+1)
            if v[i-1] == w[j-1]
                ms = s[i-1,j-1] + match_penalty
            else
                ms = s[i-1,j-1] + mismatch_penalty
            end
            test = [ms, s[i-1,j] + deletion_penalty, s[i,j-1] + deletion_penalty]
            p = argmax(test)
            s[i,j] = test[p]
            b[i,j] = p
        end
    end

    i = n1+1
    j = n2+1
    sv = []
    sw = []
    while(i>1 || j>1)
        p = b[i,j]
        if (p == 1)
            i = i-1
            j = j-1
            push!(sv, v[i])
            push!(sw, w[j])
        elseif p == 2
            i = i-1
            push!(sv, v[i])
            push!(sw, "-")
        elseif p == 3
            j = j-1
            push!(sv, "-")
            push!(sw, w[j])
        else
            break
        end
    end

    return (s[n1+1, n2+1], join(reverse(sv)), join(reverse(sw)))
end

#@time global_align("ACGTAAAAAA","GGNNNNATATATATAT")