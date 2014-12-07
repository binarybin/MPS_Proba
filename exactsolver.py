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

    def Interpreter(self):
    	if type(self.model) != AngryBoys:
    		raise Exception("The model is not supported")

    	#read vars- intitial state, state size
    	state_size = self.model.size 
    	p  = self.model.remain_proba 

    	#create sigma_x matrix
    	sigmax = np.matrix('0 1; 1 0')

    	#create hamiltonian
    	H = p*np.identity(math.pow(2,state_size))
    	matrix = 1
    	part = np.zeros((math.pow(2,state_size),math.pow(2,state_size)))
    	for i in range(1,state_size):
	    	for j in range(1,state_size):
	    		if j != i and j != i+1:
		    		matrix = (1-p)/(state_size-1)*np.kron(matrix, np.identity(2))
		    	else:
			    	matrix = (1-p)/(state_size-1)*np.kron(matrix, sigmax)
			
			part = np.add(part, matrix)
			matrix = 1
        H = np.add(H,part)

    	#call step
    	Step(H)

    def Step(self,H):
    	#take step
    	new_state = H.dot(self.model.init_state) 

    	#call evolve
    	Evolve(new_state)

    def Evolve(self, new_state):
    	#add state to list
    	self.output1.append(new_state)

    	#aux data with state (sum to 1)
    	check_norm = new_state.sum(axis=0)
    	self.output2.append(check_norm)
