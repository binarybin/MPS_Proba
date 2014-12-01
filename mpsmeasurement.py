"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpsmeasurement.py
File purpose: the class for measuring several quantities using the exact solution
Responsible person: Peiqi Wang
"""

from measurement import Measurement
from exactsolver import ExactSolver

class MpsMeasurement(Measurement):
    """
    The Measurement class for the mps method
    """
    def measureProba(self, task):
        """
        This implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
        """
        raise NotImplementedError("Please implement")
        
    def measureCorrelation(self, task):
        """
        This implements the measurement of the n-point correlation function
        """
        raise NotImplementedError("Please implement")
        
    def measureMean(self, task):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        raise NotImplementedError("please implement")
        
    def measureVariance(self, task):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        raise NotImplementedError("please implement")
    