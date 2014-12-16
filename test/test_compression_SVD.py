#compression parameter d : maximal dimension for auxillary space
import numpy as np
import copy
from cmpSVD_for_test import MpsSolver 

d=3;

L=5;
mps=[];
new_mps=np.zeros(shape=(2,1,4),dtype=float);
mps.append(new_mps)
mps[0][0,0,3]=1
mps[0][1,0,2]=1
for i in range(2,L):
    #new_mps=np.zeros(shape=(2,4,4),dtype=float);
    new_mps=np.random.rand(2,4,4);
    mps.append(new_mps)
    mps[i-1][0,1,1]=1;
    mps[i-1][0,2,3]=2;
    mps[i-1][1,3,1]=1;
    mps[i-1][1,3,3]=5;

#new_mps=np.zeros(shape=(2,4,1),dtype=float)
new_mps=np.random.rand(2,4,1);
mps.append(new_mps)
mps[L-1][0,3,0]=4
mps[L-1][1,2,0]=1


#calculate module

MpsTest=MpsSolver(mps, d)
MpsTest.CompressionSVDSweepToRight()
MpsTest.LeftNormalize(4)
print mps[1]
print '\n'
print MpsTest.mpsc[1]
print '\n'


MpsTest.CompressionSVDSweepToLeft()
print MpsTest.mpsc[1]