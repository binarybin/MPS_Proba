#compression parameter d : maximal dimension for auxillary space
import numpy as np

d=3;

L=10;
mps=[];
#new_mps=np.zeros(shape=(2,1,4),dtype=float);
new_mps=np.random.rand(2,1,4);
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

mps_module=np.tensordot(mps[0],mps[0],axes=([0],[0]))
mps_module=mps_module[0,:,0,:]
for i in range(L-1):
    B=np.tensordot(mps[i+1],mps[i+1],axes=([0],[0]))
    mps_module=np.tensordot(mps_module,B,axes=([0,1],[0,2]))

#candidate compressed MPS
mpsc=[];
#new_mps=np.zeros(shape=(2,1,d),dtype=float);
new_mps=np.random.rand(2,1,d);
mpsc.append(new_mps)
mpsc[0][0,0,2]=1
mpsc[0][1,0,0]=1
for i in range(2,L):
    #new_mps=np.zeros(shape=(2,d,d),dtype=float);
    new_mps=np.random.rand(2,d,d);
    mpsc.append(new_mps)
    mpsc[i-1][0,1,1]=1;
    mpsc[i-1][0,2,0]=2;
    mpsc[i-1][1,2,1]=1;

#new_mps=np.zeros(shape=(2,d,1),dtype=float)
new_mps=np.random.rand(2,d,1);
mpsc.append(new_mps)
mpsc[L-1][0,2,0]=1
mpsc[L-1][1,0,0]=1

##right-normalize the state
l=0;
for i in range(L-l-1):
    A=np.swapaxes(mpsc[L-1-i],1,2)
    A=np.reshape(A,(A.shape[0]*A.shape[1],A.shape[2]))
    U, s, V=np.linalg.svd(A, full_matrices=False)
    mpsc[L-1-i]=np.reshape(U,(mpsc[L-1-i].shape[0],mpsc[L-1-i].shape[2],U.shape[1]))
    mpsc[L-1-i]=np.swapaxes(mpsc[L-1-i],1,2)
    B=np.dot(np.diag(s),V)
    mpsc[L-2-i]=np.tensordot(mpsc[L-2-i],B,axes=([2],[1]))

#Initialize partial overlap list of mps and its compressed version
#In partial_overlap_lr we store the list [,[[,[[[,...., with the last element as the full overlap
#In partial_overlap_rl we store the list ,....,]]],]],] with the first element as the full overlap
partial_overlap_lr=[None]*L;
partial_overlap_rl=[None]*L;
partial_overlap_lr[0]=np.tensordot(mpsc[0],mps[0],axes=([0],[0]));
partial_overlap_lr[0]=partial_overlap_lr[0][0,:,0,:];
partial_overlap_rl[L-1]=np.tensordot(mpsc[L-1],mps[L-1],axes=([0],[0]));
partial_overlap_rl[L-1]=partial_overlap_rl[L-1][:,0,:,0];
for i in range(L-1):
    A=np.tensordot(mpsc[i+1],mps[i+1],axes=([0],[0]))
    partial_overlap_lr[i+1]=np.tensordot(partial_overlap_lr[i],A,axes=([0,1],[0,2]))

for i in range(L-1):
    A=np.tensordot(mpsc[L-2-i],mps[L-2-i],axes=([0],[0]))
    partial_overlap_rl[L-2-i]=np.tensordot(A,partial_overlap_rl[L-1-i],axes=([1,3],[0,1]))


#Perform a single sweep from left to right, return the error at the end
#CprVarSweepLR
A=np.tensordot(mps[0],partial_overlap_rl[1],axes=([2],[1]))
#perform left normalization
A=np.reshape(A,(A.shape[0],A.shape[2]))
U, s, V=np.linalg.svd(A, full_matrices=False)
mpsc[0]=np.reshape(U,(A.shape[0],1,U.shape[1]))
##Update partialoverlap list (direction left to right)
partial_overlap_lr[0]=np.tensordot(mpsc[0],mps[0],axes=([0],[0]))
partial_overlap_lr[0]=partial_overlap_lr[0][0,:,0,:]
mpsc_module=np.tensordot(mpsc[0],mpsc[0],axes=([0],[0]))
mpsc_module=mpsc_module[0,:,0,:]

for i in range(L-2):
    A=np.tensordot(mps[i+1],partial_overlap_rl[i+2],axes=([2],[1]))
    A=np.tensordot(partial_overlap_lr[i],A,axes=([1],[1]))
    A=np.swapaxes(A,0,1)
    #perform left normalization
    mA=A.shape[0]
    nA=A.shape[1]
    A=np.reshape(A,(mA*nA,A.shape[2]))
    U, s, V=np.linalg.svd(A, full_matrices=False)
    mpsc[i+1]=np.reshape(U,(mA,nA,U.shape[1]))
    ##Update partialoverlap list (direction left to right)
    A=np.tensordot(mpsc[i+1],mps[i+1],axes=([0],[0]))
    partial_overlap_lr[i+1]=np.tensordot(partial_overlap_lr[i],A,axes=([0,1],[0,2]))
    B=np.tensordot(mpsc[i+1],mpsc[i+1],axes=([0],[0]))
    mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[0,2]))

A=np.tensordot(partial_overlap_lr[L-2],mps[L-1],axes=([1],[1]))
mpsc[L-1]=np.swapaxes(A,0,1)
#no need to left-normalize the right-most MPS, we update partialoverlap list
#this gives us <mpsc,mps> at the end of the iteration
A=np.tensordot(mpsc[L-1],mps[L-1],axes=([0],[0]))
partial_overlap_lr[L-1]=np.tensordot(partial_overlap_lr[L-2],A,axes=([0,1],[0,2]))
B=np.tensordot(mpsc[L-1],mpsc[L-1],axes=([0],[0]))
mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[0,2]))
print(mps_module+mpsc_module-2*partial_overlap_lr[L-1])


#Perform a single sweep from right to left, return the error at the end
#CprVarSweepRL
A=np.tensordot(mps[L-1],partial_overlap_lr[L-2],axes=([1],[1]))
#perform right normalization
mA=A.shape[0]
nA=A.shape[1]
A=np.reshape(A,(mA*nA,A.shape[2]))
U, s, V=np.linalg.svd(A, full_matrices=False)
U=np.reshape(U,(mA,nA,U.shape[1]))
mpsc[L-1]=np.swapaxes(U,1,2)
##Update partialoverlap list (direction right to left)
partial_overlap_rl[L-1]=np.tensordot(mpsc[L-1],mps[L-1],axes=([0],[0]))
partial_overlap_rl[L-1]=partial_overlap_rl[L-1][:,0,:,0]
mpsc_module=np.tensordot(mpsc[L-1],mpsc[L-1],axes=([0],[0]))
mpsc_module=mpsc_module[:,0,:,0]

for i in range(L-2):
    A=np.tensordot(mps[L-2-i],partial_overlap_rl[L-1-i],axes=([2],[1]))
    A=np.tensordot(A,partial_overlap_lr[L-3-i],axes=([1],[1]))
    #perform right normalization
    mA=A.shape[0]
    nA=A.shape[1]
    A=np.reshape(A,(mA*nA,A.shape[2]))
    U, s, V=np.linalg.svd(A, full_matrices=False)
    U=np.reshape(U,(mA,nA,U.shape[1]))
    mpsc[L-i-2]=np.swapaxes(U,1,2)
    ##Update partialoverlap list (direction right to left)
    A=np.tensordot(mpsc[L-i-2],mps[L-i-2],axes=([0],[0]))
    partial_overlap_rl[L-i-2]=np.tensordot(partial_overlap_rl[L-i-1],A,axes=([0,1],[1,3]))
    B=np.tensordot(mpsc[L-i-2],mpsc[L-i-2],axes=([0],[0]))
    mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[1,3]))

mpsc[0]=np.tensordot(mps[0],partial_overlap_rl[1],axes=([2],[1]))
#no need to left-normalize the right-most MPS, we update partialoverlap list
#this gives us <mpsc,mps> at the end of the iteration
A=np.tensordot(mpsc[0],mps[0],axes=([0],[0]))
partial_overlap_rl[0]=np.tensordot(partial_overlap_rl[1],A,axes=([0,1],[1,3]))
B=np.tensordot(mpsc[0],mpsc[0],axes=([0],[0]))
mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[1,3]))
print(mps_module+mpsc_module-2*partial_overlap_rl[0])