"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: projectionboys.py
File purpose: the projection boys model
Responsible person: Bin Xu
"""

from numpy import zeros
import numpy as np
from model import Model

class ProjectionBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    """

    def __init__(self, size, p0, p1, q1, q2, init_state):
        super(ProjectionBoys, self).__init__()
        self.size = size
        self.p0 = p0
        self.p1 = p1
        self.q1 = q1
        self.q2 = q2
        self.init_state = init_state
        self.model_type = "ProjectionBoys"
        self.hamiltonian = r"H = p0 I + \sum_{i=1}^{n-1}\frac{p 1}{n-1}\sigma_i^x\otimes\sigma_{i+1}^x + \frac{q 1}{n-1}\pi_i^+\otimes\pi_{i+1}^- + \frac{q 2}{n-1}\pi_i^+\otimes\pi_{i+1}^-"
        self.normalizeTMat()
        
    def normalizeTMat(self):
        totalproba = self.p0 + self.p1 + self.q1 + self.q2
        self.p0 /= totalproba
        self.p1 /= totalproba
        self.p1 /= self.size - 1
        self.q1 /= totalproba
        self.q1 /= self.size - 1 
        self.q2 /= totalproba
        self.q2 /= self.size - 1

    def prepareMpo(self):
        #initialize the MPO
        self.mpo = []

        mpo_left = zeros(shape = (2, 2, 1, 5), dtype = float)
        mpo_middle = zeros(shape = (2, 2, 5, 5), dtype = float)
        mpo_right = zeros(shape = (2, 2, 5, 1), dtype = float)

        # remember our convention: phys_in, phys_out, aux_l, aux_r
        # mpo_left = [p0 I, p1 Sx, q1 Pi+, q2 Pi-, I]

        mpo_left[:, :, 0, 0] = self.p0 * self.I
        mpo_left[:, :, 0, 1] = self.p1 * self.sigma_x
        mpo_left[:, :, 0, 2] = self.q1 * self.pi_plus
        mpo_left[:, :, 0, 3] = self.q2 * self.pi_minus
        mpo_left[:, :, 0, 4] = self.I

        # mpo_middle = [I, 0, 0, 0, 0]
        #              [Sx, 0, 0, 0, 0]
        #              [pi+, 0, 0, 0, 0]
        #              [pi-, 0, 0, 0, 0]
        #              [0, p1 Sx, q1 pi+, q2 pi-, I]
        mpo_middle[:, :, 0, 0] = self.I
        mpo_middle[:, :, 1, 0] = self.sigma_x
        mpo_middle[:, :, 2, 0] = self.pi_plus
        mpo_middle[:, :, 3, 0] = self.pi_minus
        mpo_middle[:, :, 4, 1] = self.p1 * self.sigma_x
        mpo_middle[:, :, 4, 2] = self.q1 * self.pi_plus
        mpo_middle[:, :, 4, 3] = self.q2 * self.pi_minus
        mpo_middle[:, :, 4, 4] = self.I
        
        # mpo_right = [I, Sx, pi+, pi-, 0].transpose

        mpo_right[:, :, 0, 0] = self.I
        mpo_right[:, :, 1, 0] = self.sigma_x
        mpo_right[:, :, 2, 0] = self.pi_plus
        mpo_right[:, :, 3, 0] = self.pi_minus

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
    	pi_plus = np.matrix(self.pi_plus)
    	pi_minus = np.matrix(self.pi_minus)

    	#non changing channel
    	self.H = self.p0*np.identity(2**self.size) # not changing states

    	# sigma_x channel
        for i in range(self.size-1):
            Tmatrix = np.identity(1)
            for j in range(self.size):
                if j == i or j == i+1:
                    Tmatrix = np.kron(Tmatrix, sigmax)
                else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))
            self.H = np.add(self.H, Tmatrix * self.p1)
            
        # pi+ channel
        for i in range(self.size-1):
            Tmatrix = np.identity(1)
            for j in range(self.size):
                if j == i or j == i+1:
                    Tmatrix = np.kron(Tmatrix, pi_plus)
                else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))
            self.H = np.add(self.H, Tmatrix * self.q1)
            
        # pi- channel
        for i in range(self.size-1):
            Tmatrix = np.identity(1)
            for j in range(self.size):
                if j == i or j == i+1:
                    Tmatrix = np.kron(Tmatrix, pi_minus)
                else:
                    Tmatrix = np.kron(Tmatrix, np.identity(2))
            self.H = np.add(self.H, Tmatrix * self.q2)

    def prepareExactInitState(self):
        self.init_exact = np.zeros((2**self.size, 1))
        if self.init_state == "all down":
            self.init_exact[0] = 1
        else:
            raise Exception("Init state not supported!")

    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nremain_proba = "+str(self.remain_proba) +"\ninitial state: "+self.init_state)

