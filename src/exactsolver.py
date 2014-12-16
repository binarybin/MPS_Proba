"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: exactsolver.py
File purpose: the exact transition matrix solver class
Responsible person: Peter Gjeltema
"""
from solver import Solver
import numpy as np
from math import pow
from model import AngryBoys

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
        self.init_state = np.zeros(1)
        self.results = []        
        self.check_norm = []
        self.interpreter()

    def interpreter(self):
        if self.model.model_type == "AngryBoys":
            self.init_state = np.zeros((2**self.model.size, 1))
            if self.model.init_state == "all down":
                self.init_state[0] = 1
            else:
                raise Exception("Init state not supported!")
        
    	else:
    		raise Exception("The model is not supported")
        self.results.append(self.init_state)
    	#read vars- intitial state, state size
    	state_size = self.model.size 
    	p  = self.model.remain_proba 

    	#create sigma_x matrix
    	sigmax = np.matrix('0 1; 1 0')

    	#create hamiltonian
    	self.H = p*np.identity(2**state_size) # not changing states
    	
    	# changing states
    	matrix = np.identity(1)
    	part = np.zeros((2**state_size, 2**state_size))
    	
    	for i in range(state_size-1):
	    	for j in range(state_size):
	    		if j != i and j != i+1:
		    		matrix = np.kron(matrix, np.identity(2))
		    	else:
			    	matrix = np.kron(matrix, sigmax)
			
		part = np.add(part, matrix)
		matrix = np.identity(1)
	
	# add them
        self.H = np.add(self.H, (1-p)/(state_size-1)*part)

    def step(self):
        self.t += 1
    	return self.H.dot(self.results[-1])

    def evolve(self, nstep):
    	#add state to list
        for i in range(nstep):
            new_state = self.step()
            self.results.append(new_state)
            self.check_norm.append(new_state.sum(axis=0))