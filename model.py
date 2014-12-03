"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: model.py
File purpose: the model definition
Responsible person: Bin Xu
"""

from numpy import ndarray, zeros

class Model(object):
    """
    The abstract model class
    """


class AngryBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    """
    
    def __init__(self, size, remain_proba, init_state, output1, output2):
        self.size = size
        self.remain_proba = remain_proba
        self.init_state = init_state
        self.output1 = output1
        self.output2 = output2
        self.hamiltonian = r"H = pI + (1-p)\sum_{i=1}^{n-1}\frac{1}{n-1}\sigma_i^x\otimes\sigma_{i+1}^x"
        self.prepareMpo()
        self.prepareMps()

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
            
        mpo_left[0, 0, 0, 0] = p
        mpo_left[1, 1, 0, 0] = p
        mpo_left[1, 0, 0, 1] = q
        mpo_left[0, 1, 0, 1] = q
        mpo_left[1, 1, 0, 2] = 1
        mpo_left[0, 0, 0, 2] = 1
            
        # mpo_middle = [I, 0, 0]
        #              [Sx, 0, 0]
        #              [0, qSx, I]
        mpo_middle[1, 0, 2, 1] = q
        mpo_middle[0, 1, 2, 1] = q
        mpo_middle[1, 1, 2, 2] = 1
        mpo_middle[0, 0, 2, 2] = 1
        mpo_middle[0, 0, 0, 0] = 1
        mpo_middle[1, 1, 0, 0] = 1
        mpo_middle[1, 0, 1, 0] = 1
        mpo_middle[0, 1, 1, 0] = 1
            
        # mpo_right = [I, Sx, 0].transpose
            
        mpo_right[0, 0, 0, 0] = 1
        mpo_right[1, 1, 0, 0] = 1
        mpo_right[1, 0, 1, 0] = 1
        mpo_right[0, 1, 1, 0] = 1
            
        # store the list of mpo's
            
        self.mpo.append(mpo_left)
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
        
    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nt = "+str(self.remain_proba) +"\ninitial state: "+self.init_state)
    
    