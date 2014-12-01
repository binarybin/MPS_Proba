"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: model.py
File purpose: the model definition
"""

class Model(object):
    """
    The abstract model class
    """


class AngryBoys(Model):
    """
    A probablistic model that describes the model in human language and gives some parameters.
    """
    
    def __init__(self, size, change_rate, init_state, output1, output2):
        self.size = size
        self.change_rate = change_rate
        self.init_state = init_state
        self.output1 = output1
        self.output2 = output2
        self.hamiltonian = "Angry Boy"
        
    def __repr__(self):
        return ( "Hamiltonian: "+self.hamiltonian + "\nSystem length = "+str(self.size)+"\nt = "+str(self.change_rate) +"\ninitial state: "+self.init_state)
    
    