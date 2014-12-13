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
        self.results = []
        self.results.append(self.model.init_state)
        self.check_norm = []


    def Interpreter(self):
    	if type(self.model) != AngryBoys:
    		raise Exception("The model is not supported")

    	#read vars- intitial state, state size
    	state_size = self.model.size 
    	p  = self.model.remain_proba 

    	#create sigma_x matrix
    	sigmax = np.matrix('0 1; 1 0')

    	#create hamiltonian
    	self.H = p*np.identity(math.pow(2,state_size))
    	matrix = 1
    	part = np.zeros((math.pow(2,state_size),math.pow(2,state_size)))
    	for i in range(1,state_size):
	    	for j in range(1,state_size):
	    		if j != i and j != i+1:
		    		matrix = np.kron(matrix, np.identity(2))
		    	else:
			    	matrix = np.kron(matrix, sigmax)
			
			part = np.add(part, matrix)
			matrix = 1
        self.H = (1-p)/(state_size-1)*np.add(H,part)

    def Step(self):
        self.t += 1
    	return self.H.dot(self.results[-1])

    def Evolve(self, nstep):
    	#add state to list
        for i in range(num_steps):
            new_state = self.Step()
            self.results.append(new_state)
            self.check_norm.append(new_state.sum(axis=0))