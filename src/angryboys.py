"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: angryboys.py
File purpose: the angry boys model
Responsible person: Bin Xu
"""

from numpy import ndarray, zeros
import numpy as np
from model import Model

class AngryBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    """

    def __init__(self, size, remain_proba, init_state):
        super(AngryBoys, self).__init__()
        self.size = size
        self.remain_proba = remain_proba
        self.init_state = init_state
        self.model_type = "AngryBoys"
        self.hamiltonian = r"H = pI + (1-p)\sum_{i=1}^{n-1}\frac{1}{n-1}\sigma_i^x\otimes\sigma_{i+1}^x"


    def prepareMpo(self):
        #initialize the MPO
        self.mpo = []

        mpo_left = zeros(shape = (2, 2, 1, 3), dtype = float)
        mpo_middle = zeros(shape = (2, 2, 3, 3), dtype = float)
        mpo_right = zeros(shape = (2, 2, 3, 1), dtype = float)

        q = (1.0-self.remain_proba)/(self.size - 1.0)
        p = self.remain_proba

        # remember our convention: phys_in, phys_out, aux_l, aux_r
        # mpo_left = [pI, qSx, I]

        mpo_left[:, :, 0, 0] = p * self.I
        mpo_left[:, :, 0, 1] = q * self.sigma_x
        mpo_left[:, :, 0, 2] = self.I

        # mpo_middle = [I, 0, 0]
        #              [Sx, 0, 0]
        #              [0, qSx, I]
        mpo_middle[:, :, 2, 1] = q * self.sigma_x
        mpo_middle[:, :, 2, 2] = self.I
        mpo_middle[:, :, 0, 0] = self.I
        mpo_middle[:, :, 1, 0] = self.sigma_x

        # mpo_right = [I, Sx, 0].transpose

        mpo_right[:, :, 0, 0] = self.I
        mpo_right[:, :, 1, 0] = self.sigma_x

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
    	sigmax = np.matrix(self.sigma_x)

    	#non changing channel
    	self.H = self.remain_proba*np.identity(2**self.size) # not changing states

    	# nearest-neighbour changing channel
        for i in range(self.size-1):
            Tmatrix = np.identity(1)
            for j in range(self.size):
                if j == i or j == i+1:
                    Tmatrix = np.kron(Tmatrix, sigmax)
                else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))
            self.H = np.add(self.H, Tmatrix * (1-self.remain_proba)/(self.size - 1))

#    	for i in range(self.size-1):
#    	    Tmatrix = np.identity(1)
#	    for j in range(self.size):
#	    	  if j == i or j == i+1:
#                    Tmatrix = np.kron(Tmatrix, sigmax)
#	    	  else:
#                    Tmatrix = np.kron(Tmatrix, np.identity(2))
#        self.H = np.add(self.H, Tmatrix * (1-self.remain_proba)/(self.size -1))
#        print self.H

    def prepareExactInitState(self):
        self.init_exact = np.zeros((2**self.size, 1))
        if self.init_state == "all down":
            self.init_exact[0] = 1
        else:
            raise Exception("Init state not supported!")

    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nremain_proba = "+str(self.remain_proba) +"\ninitial state: "+self.init_state)

