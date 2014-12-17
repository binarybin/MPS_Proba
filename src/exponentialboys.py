"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: exponentialboys.py
File purpose: the exponential boys model
Responsible person: Bin Xu
"""

from numpy import ndarray, zeros
import numpy as np
from model import Model

class ExponentialBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    The model is similar to AngryBoys model but with exponentially decaying interaction.
    """
    
    def __init__(self, size, J, K, init_state, output1, output2):
        self.size = size
        self.init_state = init_state
        self.output1 = output1
        self.output2 = output2
        self.K = K
        self.J = J
        self.model_type = "ExponentialBoys"
        self.hamiltonian = r"H = P (I + J \sum_{i=1}^{n-1} \sum_{j=i+1}^{n-1} K^{j-i}\sigma_i^x\otimes\sigma_j^x)"
        self.normalizeTMat()
        self.prepareMpo()
        self.prepareMps()
        self.prepareTransitionalMat()
        self.prepareExactInitState()
    
    def normalizeTMat(self):
        totalproba = 1.0 # the not changing channel
        for i in range(self.size):
            for j in range(i+1, self.size):
                totalproba += self.J * self.K ** (j-i)
        self.p0 = 1.0/totalproba
        self.J /= totalproba
        
    def prepareMpo(self):
        #initialize the MPO
        self.mpo = []
        
        mpo_left = zeros(shape = (2, 2, 1, 3), dtype = float)            
        mpo_middle = zeros(shape = (2, 2, 3, 3), dtype = float)
        mpo_right = zeros(shape = (2, 2, 3, 1), dtype = float)
            
        # remember our convention: phys_in, phys_out, aux_l, aux_r
        # mpo_left = [p0 I, J K Sx, I]   
                    
        mpo_left[0, 0, 0, 0] = self.p0
        mpo_left[1, 1, 0, 0] = self.p0
        mpo_left[1, 0, 0, 1] = self.J * self.K
        mpo_left[0, 1, 0, 1] = self.J * self.K
        mpo_left[1, 1, 0, 2] = 1
        mpo_left[0, 0, 0, 2] = 1
            
        # mpo_middle = [I,  0,      0]
        #              [Sx, K I,    0]
        #              [0,  J K Sx, I]
        mpo_middle[0, 0, 0, 0] = 1
        mpo_middle[1, 1, 0, 0] = 1
        mpo_middle[1, 0, 1, 0] = 1
        mpo_middle[0, 1, 1, 0] = 1
        mpo_middle[0, 0, 1, 1] = self.K
        mpo_middle[1, 1, 1, 1] = self.K
        mpo_middle[1, 0, 2, 1] = self.J * self.K
        mpo_middle[0, 1, 2, 1] = self.J * self.K
        mpo_middle[1, 1, 2, 2] = 1
        mpo_middle[0, 0, 2, 2] = 1
            
        # mpo_right = [I, Sx, 0].transpose
            
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
	
	for i in range(self.size):  #site 1
            for j in range(i+1, self.size): #site 2, rhs of site 1
                Tmatrix = np.identity(1)
                for k in range(self.size):
                    if k == i or k == j:
                        Tmatrix = np.kron(Tmatrix, sigmax)    	
                    else:
                        Tmatrix = np.kron(Tmatrix, np.identity(2))
                self.H = np.add(self.H, Tmatrix * self.J * self.K ** (j-i))
                        
    def prepareExactInitState(self):
        self.init_exact = np.zeros((2**self.size, 1))
        if self.init_state == "all down":
            self.init_exact[0] = 1
        else:
            raise Exception("Init state not supported!")
        
    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nremain_proba = "+str(self.remain_proba) +"\ninitial state: "+self.init_state)
    
    