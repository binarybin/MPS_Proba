"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: solver.py
File purpose: the abstract solver class
"""

class Solver(object):
    """
    The abstract solver base class
    """
    
    def __init__(self, model, output1, output2):
        self.model = model
        self.output1 = output1
        self.output2 = output2        
    def Interpreter(self):
        raise NotImplementedError("Unimplemented abstract method")
    def Step(self):
        raise NotImplementedError("Unimplemented abstract method")
    def Evolve(self):
        raise NotImplementedError("Unimplemented abstract method")