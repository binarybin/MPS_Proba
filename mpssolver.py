"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpssolver.py
File purpose: the solver class based on the matrix product states
Responsible persons: 
    Liangsheng Zhang and Jun Xiong for Contraction and Compression
    Bin Xu for Interpreter
"""
from solver import Solver

class MpsSolver(Solver):
    """
    The exact solver base class, using transition matrix method
    """
          
    def Interpreter(self):
        raise NotImplementedError("please implement")
    def Step(self):
        raise NotImplementedError("please implement")
    def Evolve(self):
        raise NotImplementedError("please implement")
        
    def Contraction(self):
        raise NotImplementedError("please implement")
    def Compression(self):
        raise NotImplementedError("please implement")