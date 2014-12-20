"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: model.py
File purpose: the model definition
Responsible person: Bin Xu
"""

from numpy import zeros

class Model(object):
    """
    The abstract model class
    """
    def __init__(self):
        self.I = zeros(shape = (2, 2), dtype = float)
        self.sigma_x = zeros(shape = (2, 2), dtype = float)
        self.pi_plus = zeros(shape = (2, 2), dtype = float)
        self.pi_minus = zeros(shape = (2, 2), dtype = float)
        self.I[0, 0] = 1
        self.I[1, 1] = 1
        self.sigma_x[0, 1] = 1
        self.sigma_x[1, 0] = 1
        self.pi_plus[0, 0] = 1
        self.pi_plus[0, 1] = 1
        self.pi_minus[1, 0] = 1
        self.pi_minus[1, 1] = 1
    