from math import log2, floor
v = [-1,31, 8,30, -1, 7,-1,-1, 29,-1,26, 6, -1,-1, 2,-1,
	  -1,28,-1,-1, -1,19,25,-1, 5,-1,17,-1, 23,14, 1,-1,
	   9,-1,-1,-1, 27,-1, 3,-1, -1,-1,20,-1, 18,24,15,10,
	  -1,-1, 4,-1, 21,-1,16,11, -1,22,-1,12, 13,-1, 0,-1]

def f(ix):
    #since rax is used as an index, we assume `ix` is unsigned
    edx = 0 if ix == 0 else 2**(floor(log2(ix))+1)-1 
    magic = 0x4badf0d
    return v[(edx*magic & 0xffffffff) >> 26] #here's where I failed, I thought
                                             #wasn't byte-indexed, smh
def g(ix):
    #this edge case is important
    if ix == 0:
        return -1
    else:
        return 31 - floor(log2(ix))

#to make sure the pen and paper solution was correct
def pre(val):
    edx = val
    edx >>= 1
    edx = edx | val
    eax = edx
    eax >>= 2
    eax = eax | edx
    edx >>= 4
    edx = edx | eax
    eax = edx
    eax >>= 8
    eax = eax | edx
    edx = eax
    edx >>= 16
    edx = edx | eax
    return edx


#these are the valid numbers of edx
edx_vals = [0xffffffff, 0x7fffffff, 0x3fffffff, 0x1fffffff,
            0x0fffffff, 0x07ffffff, 0x03ffffff, 0x01ffffff,
            0x00ffffff, 0x007fffff, 0x003fffff, 0x001fffff,
            0x000fffff, 0x0007ffff, 0x0003ffff, 0x0001ffff,
            0x0000ffff, 0x00007fff, 0x00003fff, 0x00001fff,
            0x00000fff, 0x000007ff, 0x000003ff, 0x000001ff,
            0x000000ff, 0x0000007f, 0x0000003f, 0x0000001f,
            0x0000000f, 0x00000007, 0x00000003, 0x00000001,
            0x0]

log_vals = [31, 30, 29, 28,
            27, 26, 25, 24,
            23, 22, 21, 20,
            19, 18, 17, 16,
            15, 14, 13, 12,
            11, 10, 9, 8,
            7, 6, 5, 4,
            3, 2, 1, 0]

#the output of the below is 0, 1, 2, ..., 31, -1
for e in edx_vals:
    print(f(e))

#solution:
#ix |-> floor(log_2(ix)) |-> 31 - floor(log_2(ix)) if ix != 0 else -1
#see also g above

print(all([f(e) == g(e) for e in edx_vals]))


