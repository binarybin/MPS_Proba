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
from numpy import ndarray
from model import AngryBoys
import numpy as np

class MpsSolver(Solver):
    """
    The MPS solver base class
    """
    
    mps_element = ndarray(shape = (2, 10, 10), dtype = float) # this is just an example of the mps, the order of indices: physical, left_aux, right_aux
    mpo_element = ndarray(shape = (2, 2, 4, 4), dtype = float) # this is just an example of the mpo, the order of indices: physical_in, physical_out, left_aux, right_aux
    
    mps_result = [] # list of mps_chain, result history
    
    
    def __init__(self, model, bound_dimension, output1, output2):
        self.model = model
        self.output1 = output1
        self.output2 = output2
        self.bound_dimension = bound_dimension
        
    def Interpreter(self):
        if type(self.model) == AngryBoys:
            self.mpo = self.model.mpo
            self.mps = self.model.mps
        else:
            raise Exception("The model is not supported!")
        
        
    def Step(self):
        raise NotImplementedError("please implement")
    def Evolve(self):
        raise NotImplementedError("please implement")
            
    def Contraction(self):
        raise NotImplementedError("please implement")
    def CompressionSVD(self):
        """
        The compression based on SVD, to be implemented by Jun Xiong
        """
        raise NotImplementedError("please implement")
    def CompressionVariational(self):
        """
        The compression based on the variational principle, to be implemented by Peiqi Wang
        """
        raise NotImplementedError("please implement")
    def Compression(self):
        self.CompressionSVD()
        self.CompressionVariational()
    
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