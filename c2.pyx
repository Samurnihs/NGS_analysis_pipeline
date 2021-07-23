cpdef int comm(char x, char y): # returns 1 if variables are equal
    if x == y:
        return 1
    else:
        return 0


cpdef int co_st(char *a, char *b): # returns number of matches between 2 strings
    cdef int r = 0
    cdef int l = min(len(a), len(b))
    for i in range(min(len(a), len(b))):
        r+=comm(a[i], b[i])
    return r


cpdef list mismatch(char *a, char *b): # returns index of minimal absolute number of mismatches and absolute number of mismatches
    cdef int ind = 0 # index
    cdef int sim = 0 #  similarity = absolute number of matches
    cdef int m = 0
    for i in range(len(b) - len(a) + 1):
        m = co_st(a, b[i:])
        if m > sim:
            sim = m
            ind = i
    return [ind, sim/len(a)]

