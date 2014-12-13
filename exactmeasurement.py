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
    def __init__(self, solver, output1, output2):
        self.solver = solver
        self.basis = list(itertools.product(*[(0, 1)] * self.solver.model.size))
        self.output1 = output1
        self.output2 = output2
        self.measurement_list = []
        self.measure_result_list = []
        
    def convert(self, state):
        if state == "up":
            return 1
        elif state == "down":
            return 0
        else:
            raise Exception("Converting state failure")
    
    def measureProba(self, task):
        """
        This implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
        """
        if task[0] != "Proba":
            raise Exception("Task is wrong type")
        # ("Proba", [(1, "up"), (3, "down")])
        temptask = []
        proba = 0
        time = self.getTimeTask(task, temptask)
        for i in range(len(self.basis)):
            match = True
            for j in range(len(temptask)):
                if self.basis[i][temptask[j][0]-1] != self.convert(temptask[j][1]):
                    match = False
            if match:
                proba += self.solver.results[time][i]
        return proba

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
    