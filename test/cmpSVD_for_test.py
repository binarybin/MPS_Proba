import numpy as np
import copy

class MpsSolver:
    """
        The MPS solver base class
        """
    """
        Variables in the class
        L: length of the chain
        d: dimension of the auxilary space for the compressed MPS
        n: dimension of the physical space
        mps: current state expressed in mps. It is the mps obtained after apply the mpo
        mpsc: current compressed mps. It is obtained by calling compression algorithm
        mpo: operator. By convention we always apply mpo on mpsc
        partial_overlap_lr, partial_overlap_rl: the partial overlap between mps and mpsc. this is only useful in the algorithm of compression by variation
        mps_result: list of mps keeping the history
        t: current time
        epsil: error threshold for CompressionVariational.
        cpr_err: error of the compression(L2 distance between compressed state and true state)
        """
    ##mps_element = ndarray(shape = (2, 10, 10), dtype = float) # this is just an example of the mps, the order of indices: physical, left_aux, right_aux
    ##mpo_element = ndarray(shape = (2, 2, 4, 4), dtype = float) # this is just an example of the mpo, the order of indices: physical_in, physical_out, left_aux, right_aux
    
    mps_result = [] # list of mps_chain, result history
    
    
    def __init__(self, mps, d):
        self.mps = mps
        self.d=d
        #self.mps_fist=mps[0]
        #self.n=self.mps_fist.shape[0]
        self.n=self.mps[0].shape[0]
        self.L=len(self.mps)
    

    def CompressionSVDSweepToRight(self):
    # sweep from left to right and compress each matrix in mps
        self.mpsc=[]
        #mpsc records the result
        self.mpsc.append(self.mps[0]) # mps[0] has dimension 1*x
        mps_current=self.mps[1] #do svd and compression on site i
        for i in range(1,self.L-1):
            shape_left=mps_current.shape[1]
            shape_right=mps_current.shape[2]
            #these are the dimensions of the original matrix
            s_dim=min(self.d, min(shape_left,shape_right))
            #the dimension of the compressed s matrix
            new_mpsc=np.random.rand(self.n, shape_left, s_dim)
            mps_next=np.random.rand(self.n, s_dim, self.mps[i+1].shape[2])
            #record the altered matrix of mps[i+1]
            for j in range(0, mps_current.shape[0]):
                U, s, V=np.linalg.svd(mps_current[j], full_matrices=False)
                U=U[:, 0:s_dim]
                s1=np.diag(s)[0:s_dim, 0:s_dim]
                V=V[0:s_dim, :]
                new_mpsc[j]=U
                B=np.dot(s1,V)
                C=np.dot(B,self.mps[i+1][j])
                mps_next[j]=C
            self.mpsc.append(new_mpsc)
            mps_current=mps_next
        self.mpsc.append(mps_current)
        
    def CompressionSVDSweepToLeft(self):
        # sweep from Right to left and compress each matrix in mps
        self.mpsc=[]
        #mpsc records the result
        self.mpsc.insert(0,self.mps[self.L-1]) # mps[0] has dimension 1*x
        mps_current=self.mps[self.L-2] #do svd and compression on site i
        for i in range(0,self.L-2):
            shape_left=mps_current.shape[1]
            shape_right=mps_current.shape[2]
            #these are the dimensions of the original matrix
            s_dim=min(self.d, min(shape_left,shape_right))
            #the dimension of the compressed s matrix
            new_mpsc=np.random.rand(self.n, s_dim, shape_right)
            mps_next=np.random.rand(self.n, self.mps[self.L-3-i].shape[1], s_dim)
            #record the altered matrix of mps[i+1]
            for j in range(0, mps_current.shape[0]):
                U, s, V=np.linalg.svd(mps_current[j], full_matrices=False)
                U=U[:, 0:s_dim]
                s1=np.diag(s)[0:s_dim, 0:s_dim]
                V=V[0:s_dim, :]
                new_mpsc[j]=V
                B=np.dot(U, s1)
                C=np.dot(self.mps[self.L-3-i][j], B)
                mps_next[j]=C
            self.mpsc.insert(0,new_mpsc)
            mps_current=mps_next
        self.mpsc.insert(0,mps_current)

    def LeftNormalize(self,l):
        for i in range(0,l-1):
            A=np.reshape(self.mps[i],(self.mps[i].shape[0]*self.mps[i].shape[1],self.mps[i].shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            self.mps[i]=np.reshape(U,(self.mps[i].shape[0],self.mps[i].shape[1],U.shape[1]))
            B=np.dot(np.diag(s),V)
            self.mps[i+1]=np.tensordot(self.mps[i+1],B,axes=([1],[1]))
            self.mps[i+1]=np.swapaxes(self.mps[i+1],1,2)


