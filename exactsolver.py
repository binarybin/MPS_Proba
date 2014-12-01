"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: exactsolver.py
File purpose: the exact transition matrix solver class
"""
from solver import Solver

class ExactSolver(Solver):
    """
    The exact solver base class, using transition matrix method
    """
          
    def Interpreter(self):
        raise NotImplementedError("please implement")
    def Step(self):
        raise NotImplementedError("please implement")
    def Evolve(self):
        raise NotImplementedError("please implement")