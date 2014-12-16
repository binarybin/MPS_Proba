"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: exactsolver.py
File purpose: the exact transition matrix solver class
Responsible person: Peter Gjeltema
"""
from solver import Solver
import numpy as np
from angryboys import AngryBoys
from radiatingboys import RadiatingBoys
from copy import deepcopy

class ExactSolver(Solver):
    """
    The exact solver base class, using transition matrix method
    """
    def __init__(self,model,output1,output2):
    	self.model = model      
    	self.output1 = output1
    	self.output2 = output2
        self.H = 0
        self.t = 0
        self.results = []        
        self.check_norm = []
        self.interpreter()
    
    def interpreter(self):
        if self.model.model_type == "AngryBoys" or self.model.model_type == "RadiatingBoys":
            self.H = self.model.H
            self.state = deepcopy(self.model.init_exact)
            self.results.append(self.state)
        else:
            raise Exception("The model is not supported!")


    def step(self):
        self.t += 1
    	return self.H.dot(self.results[-1])

    def evolve(self, nstep):
    	#add state to list
        for i in range(nstep):
            new_state = self.step()
            self.results.append(new_state)
            self.check_norm.append(new_state.sum(axis=0))