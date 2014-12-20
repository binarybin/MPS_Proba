"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: radiatingboys.py
File purpose: the radiating boys model
Responsible person: Bin Xu
"""

from numpy import ndarray, zeros
import numpy as np
from model import Model

class RadiatingBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    The model is similar to AngryBoys model but with second nearest neighbour interaction.
    """
    
    def __init__(self, size, remain_proba, nearest_neighbour_proba, second_neighbour_proba, init_state):
        self.size = size
        self.p0 = remain_proba
        self.p1 = nearest_neighbour_proba
        self.p2 = second_neighbour_proba
        if self.p0 + self.p1 + self.p2 != 1:
            raise Exception("The 3 probabilities should add to 1!")
        self.init_state = init_state
        self.model_type = "RadiatingBoys"
        self.hamiltonian = r"H = p0 I + p1/(n-1) \sum_{i=1}^{n-1}\sigma_i^x\otimes\sigma_{i+1}^x + p2/(n-2) \sum_{i=1}^{n-2} \sigma_i^x\otimes\sigma_{i+2}^x"
        
    def prepareMpo(self):
        #initialize the MPO
        self.mpo = []
        
        mpo_left = zeros(shape = (2, 2, 1, 4), dtype = float)            
        mpo_middle = zeros(shape = (2, 2, 4, 4), dtype = float)
        mpo_right = zeros(shape = (2, 2, 4, 1), dtype = float)
            
        # remember our convention: phys_in, phys_out, aux_l, aux_r
        # mpo_left = [p0 I, p1 Sx, p2 Sx, I]
        
        self.p1 /= self.size-1
        self.p2 /= self.size-2    
                    
        mpo_left[0, 0, 0, 0] = self.p0
        mpo_left[1, 1, 0, 0] = self.p0
        mpo_left[1, 0, 0, 1] = self.p1
        mpo_left[0, 1, 0, 1] = self.p1
        mpo_left[1, 0, 0, 2] = self.p2
        mpo_left[0, 1, 0, 2] = self.p2
        mpo_left[1, 1, 0, 3] = 1
        mpo_left[0, 0, 0, 3] = 1
            
        # mpo_middle = [I, 0, 0, 0]
        #              [Sx, 0, 0, 0]
        #              [0, I, 0, 0]
        #              [0, p1 Sx, p2 Sx, I]
        mpo_middle[0, 0, 0, 0] = 1
        mpo_middle[1, 1, 0, 0] = 1
        mpo_middle[1, 0, 1, 0] = 1
        mpo_middle[0, 1, 1, 0] = 1
        mpo_middle[0, 0, 2, 1] = 1
        mpo_middle[1, 1, 2, 1] = 1
        mpo_middle[1, 0, 3, 1] = self.p1
        mpo_middle[0, 1, 3, 1] = self.p1
        mpo_middle[1, 0, 3, 2] = self.p2
        mpo_middle[0, 1, 3, 2] = self.p2
        mpo_middle[1, 1, 3, 3] = 1
        mpo_middle[0, 0, 3, 3] = 1
            
        # mpo_right = [I, Sx, 0, 0].transpose
            
        mpo_right[0, 0, 0, 0] = 1
        mpo_right[1, 1, 0, 0] = 1
        mpo_right[1, 0, 1, 0] = 1
        mpo_right[0, 1, 1, 0] = 1
            
        # store the list of mpo's
            
        self.mpo.append(mpo_left)
        for i in range(self.size-2):
            self.mpo.append(mpo_middle)
        self.mpo.append(mpo_right)
        
    def prepareMps(self):
        self.mps = []
        if self.init_state == "all down":
            for i in range(self.size):
                new_mps = zeros(shape = (2, 1, 1), dtype = float)
                new_mps[0, 0, 0] = 1
                self.mps.append(new_mps)
        elif type(self.init_state) == list:
            if len(self.init_state) != self.size:
                raise Exception("The size of the initial condition does not match with the size of the model.")
            for i in range(self.size):
                new_mps = zeros(shape = (2, 1, 1), dtype = float)
                if self.init_state[i] == 0:
                    new_mps[0, 0, 0] = 1
                elif self.init_state[i] == 1:
                    new_mps[1, 0, 0] = 1
                else:
                    raise Exception("Initial condition can only have 0 or 1 for this model.")
                self.mps.append(new_mps)
        else:
            raise Exception("Initial condition not supported!")
            
    def prepareTransitionalMat(self):
    	#create sigma_x matrix
    	sigmax = np.matrix('0 1; 1 0')

    	#non changing channel
    	self.H = self.p0*np.identity(2**self.size) # not changing states
    	
    	# nearest-neighbour changing channel	
    	for i in range(self.size-1):
    	    Tmatrix = np.identity(1)
	    for j in range(self.size):
	    	if j == i or j == i+1:
                    Tmatrix = np.kron(Tmatrix, sigmax)
	    	else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))	
            self.H = np.add(self.H, Tmatrix * self.p1)
		
	# second-neighbour changing channel	
    	for i in range(self.size-2):
    	    Tmatrix = np.identity(1)
	    for j in range(self.size):
	    	if j == i or j == i+2:
                    Tmatrix = np.kron(Tmatrix, sigmax)
	    	else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))	
	    self.H = np.add(self.H, Tmatrix * self.p2)
	
        
    def prepareExactInitState(self):
        self.init_exact = np.zeros((2**self.size, 1))
        if self.init_state == "all down":
            self.init_exact[0] = 1
        else:
            raise Exception("Init state not supported!")
        
    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nremain_proba = "+str(self.remain_proba) +"\ninitial state: "+self.init_state)
    
    
