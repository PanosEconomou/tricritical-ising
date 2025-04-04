import numpy as np
from fractions import Fraction

np.set_printoptions(formatter={'all': lambda x: str(Fraction(x).limit_denominator())})

BLD = '\033[1m'
RST = '\033[0m'
BRD = '\033[91m'
BCY = '\033[96m'
BGR = '\033[32m'

k   = 8;
L   = np.arange(k+1)/2
LM  = np.array([(l,m/2) for l in L for m in range(-int(2*l),int(2*l)+1,2)])
i   = lambda l,m: int(2*l * (l+1) + m)
idx = [True]*len(LM)
for (j,I) in enumerate(idx):
    if I:
        (l,m) = LM[j]
        if - k/2 + l <= m+k/2 <= k/2 - l:
            idx[i(k/2 - l,m + k/2)] = False
        if - k/2 + l <= m-k/2 <= k/2 - l:
            idx[i(k/2 - l,m - k/2)] = False

LM = LM[idx]

H   = lambda l,m,k: l*(l+1)/(k+2) - m**2/k
h   = np.array([H(*lm,k) for lm in LM])

hi = np.array([ i+j for i in [0,7/16,1/10,3/5,3/80,3/2] for j in [0,7/16,1/10,3/5,3/80,3/2]])


print(BLD+BGR+"Here are the conformal weights"+RST)
print(BLD+BCY+"\nCoset: "+RST,np.sort(h),len(h))
print(BLD+BRD+"\nIsing: "+RST,np.sort(hi), len(hi))
print()
[print([float(i),h[np.where((i-h)%1 == 0)[0]].tolist()]) for i in hi]


