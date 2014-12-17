"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: solver.py
File purpose: the abstract solver class
Responsible person: Bin Xu
"""

class Solver(object):
    """
    The abstract solver base class
    """
    boy_models = ["AngryBoys", "RadiatingBoys", "ExponentialBoys"]
    def __init__(self, model, output1, output2):
        self.model = model
        self.output1 = output1
        self.output2 = output2        
    def interpreter(self):
        raise NotImplementedError("Unimplemented abstract method")
    def step(self):
        raise NotImplementedError("Unimplemented abstract method")
    def evolve(self):
        raise NotImplementedError("Unimplemented abstract method")