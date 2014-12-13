"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpssolver.py
File purpose: the solver class based on the matrix product states
Responsible persons: 
    Peiqi Wang and Jun Xiong for Contraction and Compression
    Bin Xu for Interpreter
"""
from solver import Solver
##from model import AngryBoys
from copy import deepcopy
import numpy as np

class MpsSolver(Solver):
    """
    The MPS solver base class
    """
    """
    Variables in the class
    L: length of the chain
    bound_dimension: dimension of the auxilary space for the compressed MPS
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
    """CompressionSVDSweepToRight and CompressionSVDSweepToLeft are the two methods for SVD compression implemented by Jun Xiong. The result is stored in self.mpsc.
        """
    ##mps_element = ndarray(shape = (2, 10, 10), dtype = float) # this is just an example of the mps, the order of indices: physical, left_aux, right_aux
    ##mpo_element = ndarray(shape = (2, 2, 4, 4), dtype = float) # this is just an example of the mpo, the order of indices: physical_in, physical_out, left_aux, right_aux
    
    
    def __init__(self, model, bound_dimension, output1, output2):
        self.model = model
        self.output1 = output1
        self.output2 = output2
        self.bound_dimension = bound_dimension 
        self.t=0
        self.mps_result = [] # list of mps_chain, result history
        self.Interpreter()
        
    def Interpreter(self):
        if self.model.model_type == "AngryBoys":
            self.mpo = self.model.mpo
            self.mps = self.model.mps
            #when initializing, put mpsc the same as mps so we can apply mpo on it
            self.mpsc = deepcopy(self.model.mps)
            self.mps_result.append(self.mpsc)
            self.L = len(self.model.mps)
            self.n = np.shape(self.model.mps[0])[0]
            self.partial_overlap_lr=[None]*self.L;
            self.partial_overlap_rl=[None]*self.L;
            self.cpr_err=0
            self.epsil=0.01
        else:
            raise Exception("The model is not supported!")
            
    """def CompressionSVD(self):
        
        The compression based on SVD, to be implemented by Jun Xiong
        
        raise NotImplementedError("please implement")"""
        
    def Compression(self):
        self.CompressionSVD()
        self.CompressionVariational()

    def Evolve(self,nstep):
        for i in range(nstep):
            self.Step()
          
    #Update the system state from t to t+1
    def Step(self):
        self.t=self.t+1
        self.Contraction()
        self.CompressionVariational()
        self.mps_result.append(self.mpsc)
    
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
            s_dim=min(self.bound_dimension, min(shape_left,shape_right))
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
            s_dim=min(self.bound_dimension, min(shape_left,shape_right))
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
        
    #apply mpo on the current compressed mps (mpsc). store the result on variable mps
    #convention for mpo: phys_in, phys_out, aux_l, aux_r
    #convention for mps: phys, aux_l, aux_r
    def Contraction(self):
        for i in range(0,self.L):
            A=np.tensordot(self.mpo[i],self.mpsc[i],axes=([0],[0]))
            A=np.swapaxes(A,2,3)
            self.mps[i]=np.reshape(A,(A.shape[0], A.shape[1]*A.shape[2], A.shape[3]*A.shape[4]))
    
    #overlap two mps, output <mps1,mps2>
    def Overlap(self,mps1,mps2):
        result=np.tensordot(mps1[0],mps2[0],axes=([0],[0]))
        result=result[0,:,0,:]
        L=len(mps1)
        if len(mps2)!=L:
            raise Exception("Cannot overlap two mps with different lengths")
        for i in range(L-1):
            B=np.tensordot(mps1[i+1],mps2[i+1],axes=([0],[0]))
            result=np.tensordot(result,B,axes=([0,1],[0,2]))
        return result
    
    #left-normalize the MPS from the left end to MPS[l]
    def LeftNormalize(self,l):
        for i in range(0,l-1):
            A=np.reshape(self.mps[i],(self.mps[i].shape[0]*self.mps[i].shape[1],self.mps[i].shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            self.mps[i]=np.reshape(U,(self.mps[i].shape[0],self.mps[i].shape[1],U.shape[1]))
            B=np.dot(np.diag(s),V)
            self.mps[i+1]=np.tensordot(self.mps[i+1],B,axes=([1],[1]))
            self.mps[i+1]=np.swapaxes(self.mps[i+1],1,2)

    #right-normalize the MPS from the right end to MPS[l]
    def RightNormalize(self,l):
        L=len(self.mps)
        for i in range(L-l-1):
            A=np.swapaxes(self.mps[L-1-i],1,2)
            A=np.reshape(A,(A.shape[0]*A.shape[1],A.shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            self.mps[L-1-i]=np.reshape(U,(self.mps[L-1-i].shape[0],self.mps[L-1-i].shape[2],U.shape[1]))
            self.mps[L-1-i]=np.swapaxes(self.mps[L-1-i],1,2)
            B=np.dot(np.diag(s),V)
            self.mps[L-2-i]=np.tensordot(self.mps[L-2-i],B,axes=([2],[1]))

    #obtain a mixed-canonical form centered on MPS[l]
    def MixedCanonize(self,l):
        self.LeftNormalize(l);
        self.RightNormalize(l);
    
    '''
    The following code implements the Compression by variation.
    '''
    #Form a random guess for mpsc,and right normalized it. The result serves as the starting point for the iterations
    def InitializeMpscVar(self):
        self.mpsc=[];
        #new_mps=np.zeros(shape=(n,1,d),dtype=float);
        new_mps=np.random.rand(self.n,1,self.bound_dimension);
        self.mpsc.append(new_mps)
        for i in range(2,self.L):
            #new_mps=np.zeros(shape=(n,d,d),dtype=float);
            new_mps=np.random.rand(self.n,self.bound_dimension,self.bound_dimension);
            self.mpsc.append(new_mps)       
        #new_mps=np.zeros(shape=(2,d,1),dtype=float)
        new_mps=np.random.rand(self.n,self.bound_dimension,1);
        self.mpsc.append(new_mps)
        #right-normalizae the states
        for i in range(self.L-1):
            A=np.swapaxes(self.mpsc[self.L-1-i],1,2)
            A=np.reshape(A,(A.shape[0]*A.shape[1],A.shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            self.mpsc[self.L-1-i]=np.reshape(U,(self.mpsc[self.L-1-i].shape[0],self.mpsc[self.L-1-i].shape[2],U.shape[1]))
            self.mpsc[self.L-1-i]=np.swapaxes(self.mpsc[self.L-1-i],1,2)
            B=np.dot(np.diag(s),V)
            self.mpsc[self.L-2-i]=np.tensordot(self.mpsc[self.L-2-i],B,axes=([2],[1]))

    #Initialize the two lists of partial overlap
    def InitializePartialOvl(self):
        self.partial_overlap_lr[0]=np.tensordot(self.mpsc[0],self.mps[0],axes=([0],[0]));
        self.partial_overlap_lr[0]=self.partial_overlap_lr[0][0,:,0,:];
        self.partial_overlap_rl[self.L-1]=np.tensordot(self.mpsc[self.L-1],self.mps[self.L-1],axes=([0],[0]));
        self.partial_overlap_rl[self.L-1]=self.partial_overlap_rl[self.L-1][:,0,:,0];
        for i in range(self.L-1):
            A=np.tensordot(self.mpsc[i+1],self.mps[i+1],axes=([0],[0]))
            self.partial_overlap_lr[i+1]=np.tensordot(self.partial_overlap_lr[i],A,axes=([0,1],[0,2]))
        for i in range(self.L-1):
            A=np.tensordot(self.mpsc[self.L-2-i],self.mps[self.L-2-i],axes=([0],[0]))
            self.partial_overlap_rl[self.L-2-i]=np.tensordot(A,self.partial_overlap_rl[self.L-1-i],axes=([1,3],[0,1]))

    #Perform a single sweep from left to right
    def CompressionSweepLeftRight(self):
        A=np.tensordot(self.mps[0],self.partial_overlap_rl[1],axes=([2],[1]))
        #perform left normalization
        A=np.reshape(A,(A.shape[0],A.shape[2]))
        U, s, V=np.linalg.svd(A, full_matrices=False)
        self.mpsc[0]=np.reshape(U,(A.shape[0],1,U.shape[1]))
        ##Update partialoverlap list (direction left to right)
        self.partial_overlap_lr[0]=np.tensordot(self.mpsc[0],self.mps[0],axes=([0],[0]))
        self.partial_overlap_lr[0]=self.partial_overlap_lr[0][0,:,0,:]
        mpsc_module=np.tensordot(self.mpsc[0],self.mpsc[0],axes=([0],[0]))
        mpsc_module=mpsc_module[0,:,0,:]
        for i in range(self.L-2):
            A=np.tensordot(self.mps[i+1],self.partial_overlap_rl[i+2],axes=([2],[1]))
            A=np.tensordot(self.partial_overlap_lr[i],A,axes=([1],[1]))
            A=np.swapaxes(A,0,1)
            #perform left normalization
            mA=A.shape[0]
            nA=A.shape[1]
            A=np.reshape(A,(mA*nA,A.shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            self.mpsc[i+1]=np.reshape(U,(mA,nA,U.shape[1]))
            ##Update partialoverlap list (direction left to right)
            A=np.tensordot(self.mpsc[i+1],self.mps[i+1],axes=([0],[0]))
            self.partial_overlap_lr[i+1]=np.tensordot(self.partial_overlap_lr[i],A,axes=([0,1],[0,2]))
            B=np.tensordot(self.mpsc[i+1],self.mpsc[i+1],axes=([0],[0]))
            mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[0,2]))
            print "We are sweeping to right on site", i
        A=np.tensordot(self.partial_overlap_lr[self.L-2],self.mps[self.L-1],axes=([1],[1]))
        self.mpsc[self.L-1]=np.swapaxes(A,0,1)
        #no need to left-normalize the right-most MPS, we update partialoverlap list
        #this gives us <mpsc,mps> at the end of the iteration
        A=np.tensordot(self.mpsc[self.L-1],self.mps[self.L-1],axes=([0],[0]))
        self.partial_overlap_lr[self.L-1]=np.tensordot(self.partial_overlap_lr[self.L-2],A,axes=([0,1],[0,2]))
        B=np.tensordot(self.mpsc[self.L-1],self.mpsc[self.L-1],axes=([0],[0]))
        mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[0,2]))
        self.cpr_err=mpsc_module-2*self.partial_overlap_lr[self.L-1]
    
    #Perform a single sweep from right to left
    def CompressionSweepRightLeft(self):
        A=np.tensordot(self.mps[self.L-1],self.partial_overlap_lr[self.L-2],axes=([1],[1]))
        #perform right normalization
        mA=A.shape[0]
        nA=A.shape[1]
        A=np.reshape(A,(mA*nA,A.shape[2]))
        U, s, V=np.linalg.svd(A, full_matrices=False)
        U=np.reshape(U,(mA,nA,U.shape[1]))
        self.mpsc[self.L-1]=np.swapaxes(U,1,2)
        ##Update partialoverlap list (direction right to left)
        self.partial_overlap_rl[self.L-1]=np.tensordot(self.mpsc[self.L-1],self.mps[self.L-1],axes=([0],[0]))
        self.partial_overlap_rl[self.L-1]=self.partial_overlap_rl[self.L-1][:,0,:,0]
        mpsc_module=np.tensordot(self.mpsc[self.L-1],self.mpsc[self.L-1],axes=([0],[0]))
        mpsc_module=mpsc_module[:,0,:,0]
        for i in range(self.L-2):
            A=np.tensordot(self.mps[self.L-2-i],self.partial_overlap_rl[self.L-1-i],axes=([2],[1]))
            A=np.tensordot(A,self.partial_overlap_lr[self.L-3-i],axes=([1],[1]))
            #perform right normalization
            mA=A.shape[0]
            nA=A.shape[1]
            A=np.reshape(A,(mA*nA,A.shape[2]))
            U, s, V=np.linalg.svd(A, full_matrices=False)
            U=np.reshape(U,(mA,nA,U.shape[1]))
            self.mpsc[self.L-i-2]=np.swapaxes(U,1,2)
            ##Update partialoverlap list (direction right to left)
            A=np.tensordot(self.mpsc[self.L-i-2],self.mps[self.L-i-2],axes=([0],[0]))
            self.partial_overlap_rl[self.L-i-2]=np.tensordot(self.partial_overlap_rl[self.L-i-1],A,axes=([0,1],[1,3]))
            B=np.tensordot(self.mpsc[self.L-i-2],self.mpsc[self.L-i-2],axes=([0],[0]))
            mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[1,3]))
            print "We are sweeping to left on site", i
        self.mpsc[0]=np.tensordot(self.mps[0],self.partial_overlap_rl[1],axes=([2],[1]))
        #no need to left-normalize the right-most MPS, we update partialoverlap list
        #this gives us <mpsc,mps> at the end of the iteration
        A=np.tensordot(self.mpsc[0],self.mps[0],axes=([0],[0]))
        self.partial_overlap_rl[0]=np.tensordot(self.partial_overlap_rl[1],A,axes=([0,1],[1,3]))
        B=np.tensordot(self.mpsc[0],self.mpsc[0],axes=([0],[0]))
        mpsc_module=np.tensordot(mpsc_module,B,axes=([0,1],[1,3]))
        self.cpr_err=mpsc_module-2*self.partial_overlap_rl[0]
          
    #wrap everything up
    def CompressionVariational(self):
        """
        The compression based on the variational principle, to be implemented by Peiqi Wang
        """
        ##form a initial guess
        self.InitializeMpscVar()
        ##Initialize Partial Overlap lists
        self.InitializePartialOvl()
        ##Calculate the L2 norm of the MPS to be compressed
        mps_norm=self.Overlap(self.mps,self.mps)
        error=10000
        ##Direction of previous sweep. 0 means left to right and 1 means right to left
        last_direction=1
        while (error>self.epsil):
            if (last_direction==1):
                last_cpr_err=self.cpr_err
                self.CompressionSweepLeftRight()
                self.cpr_err=self.cpr_err+mps_norm
                print(self.cpr_err)
                error=abs(last_cpr_err-self.cpr_err)
            elif (last_direction==0):
                last_cpr_err=self.cpr_err
                self.CompressionSweepRightLeft()
                self.cpr_err=self.cpr_err+mps_norm
                print(self.cpr_err)
                error=abs(last_cpr_err-self.cpr_err)