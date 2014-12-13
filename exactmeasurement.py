"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: exactmeasurement.py
File purpose: the class for measuring several quantities using the exact solution
Responsible person: Peter Gjeltema
"""

from measurement import Measurement
from exactsolver import ExactSolver
import numpy as np
import itertools
from collections import defaultdict

class ExactMeasurement(Measurement):
    """
    The Measurement class for exact solution
    """
    def __init__(self):
        self.basis = itertools.product(*[(0, 1)] * self.solver.L)

    def measureProba(self, task):
        """
        This implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
        """
        if task[0] != "Proba":
            raise Exception("Task is wrong type")
            
        state_list = {}
        counter = 0
        for item in task[1]:
            if item[1] == "up":
                state = 1
            else:
                state = 0
            for combo in basis:
                if combo(item(0)) != state:
                    break
                elif counter == 0:
                    state_list.add[combo] = 0
                else:
                    state_list[combo] += 1 
            counter+=1
        tot_count=0
        for key in state_list:
            if state_list[key] == counter:
                tot_count+=1
        return tot_count/self.solver.L

    def measureCorrelation(self, task):
        """
        This implements the measurement of the n-point correlation function
        """
        if task[0] != "Correlation":
            raise Exception("Task is wrong type")
        
    def measureMean(self, task):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Mean":
            raise Exception("Task is wrong type")
        
        
    def measureVariance(self, task):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Variance":
            raise Exception("Task is wrong type")
    