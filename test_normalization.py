import numpy as np
L=10;
mps=[];
new_mps=np.random.rand(2,1,4);
mps.append(new_mps)
mps[0][0,0,3]=1
mps[0][1,0,2]=1
for i in range(2,L+1):
    new_mps=np.random.rand(2,4,4);
    mps.append(new_mps)
    mps[i-1][0,1,1]=1;
    mps[i-1][0,2,3]=1;
    mps[i-1][1,3,1]=1;
    mps[i-1][1,3,3]=1;

#calculate state(0,0,....0) before normalization

Y=mps[0][0,:,:]  
for i in range(L-1):
    Y = np.tensordot(Y,mps[i+1][0,:,:],axes=([1],[0]))

#Mixed-canonizer MPS with center at l=4
l=4;

for i in range(l):
    A=np.reshape(mps[i],(mps[i].shape[0]*mps[i].shape[1],mps[i].shape[2]))
    U, s, V=np.linalg.svd(A, full_matrices=False)
    mps[i]=np.reshape(U,(mps[i].shape[0],mps[i].shape[1],U.shape[1]))
    B=np.dot(np.diag(s),V)
    mps[i+1]=np.tensordot(mps[i+1],B,axes=([1],[1]))
    mps[i+1]=np.swapaxes(mps[i+1],1,2)
          

for i in range(L-l-1):
    A=np.swapaxes(mps[L-1-i],1,2)
    A=np.reshape(A,(A.shape[0]*A.shape[1],A.shape[2]))
    U, s, V=np.linalg.svd(A, full_matrices=False)
    mps[L-1-i]=np.reshape(U,(mps[L-1-i].shape[0],mps[L-1-i].shape[2],U.shape[1]))
    mps[L-1-i]=np.swapaxes(mps[L-1-i],1,2)
    B=np.dot(np.diag(s),V)
    mps[L-2-i]=np.tensordot(mps[L-2-i],B,axes=([2],[1]))

#calculate state(0,0,....0) after normalization

X=mps[0][0,:,:]  
for i in range(L-1):
    X = np.tensordot(X,mps[i+1][0,:,:],axes=([1],[0]))

#test if the state(0,0,....0) before and after normalization are equal
print(np.allclose(X,Y))

#test normalization
A=np.tensordot(mps[1],mps[1],axes=([0,1],[0,1]))
B=np.tensordot(mps[L-1],mps[L-1],axes=([0,2],[0,2]))
print(A)
print(B)