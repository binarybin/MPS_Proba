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

class MpsSolver(Solver):
    """
    The exact solver base class, using transition matrix method
    """
    
    mps_element = ndarray(shape = (2, 10, 10), dtype = float) # this is just an example of the mps, the order of indices: physical, left_aux, right_aux
    mpo_element = ndarray(shape = (2, 2, 4, 4), dtype = float) # this is just an example of the mpo, the order of indices: physical_in, physical_out, left_aux, right_aux
    
    mps_chain = [] # list of mps
    
    mps_result = [] # list of mps_chain, result history
    
    
    def __init__(self, model, bound_dimension, output1, output2):
        self.model = model
        self.output1 = output1
        self.output2 = output2
        self.bound_dimension = bound_dimension
    def Interpreter(self):
        raise NotImplementedError("please implement")
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