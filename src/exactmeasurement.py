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

class ExactMeasurement(Measurement):
    """
    The Measurement class for exact solution
    """
    def __init__(self, solver):
        self.solver = solver
        self.basis = list(itertools.product(*[(0, 1)] * self.solver.model.size))
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
        #check for correct task
        if task[0] != "Proba":
            raise Exception("Task is wrong type")
        #get tasks
        temptask = []
        time = self.getTimeTask(task, temptask)

        #negative sites case
        temptask = [(self.solver.model.size + one_task[0], one_task[1]) if one_task[0] < 0 else one_task for one_task in temptask]
        
        #compute probabilities
        proba = 0
        for state_idx, basis in enumerate(self.basis):
            match = True
            for one_task in temptask:
                if basis[one_task[0]-1] != self.convert(one_task[1]):
                    match = False
            if match:
                proba += self.solver.results[time][state_idx]
        return float(proba)

    def measureCorrelation(self, task, up=None, down=None):
        """
        This implements the measurement of the n-point correlation function
        """
        if task[0] != "Correlation":
            raise Exception("Task is wrong type")

        if up is None:
            up = 1
        if down is None:
            down = -1

        #get time and tasks
        temptask = []
        time = self.getTimeTask(task, temptask)
        
        temptask2 = [self.solver.model.size + one_task if one_task < 0 else one_task for one_task in temptask] 
        temptask = temptask2        
        corr = 0
        
        for state_idx, basis in enumerate(self.basis):
            temp_corr = self.solver.results[time][state_idx]
            for one_task in temptask:
                if basis[one_task-1] == 1: temp_corr *= up # "up"
                elif basis[one_task-1] == 0: temp_corr *= down # "down"
                else: raise Exception("basis contains illegal one body states")
            corr += temp_corr
        return float(corr)
                
    def measureMean(self, task, up=None,down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Mean":
            raise Exception("Task is wrong type")

        newtask = ("Correlation", task[1], task[2])
        return self.measureCorrelation(newtask, up, down)
        
        
    def measureVariance(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Variance":
            raise Exception("Task is wrong type")

        mean_task = ("Mean", task[1], task[2])        
        average = self.measureMean(mean_task, up, down)

        square_average = self.measureMean(mean_task, 1, 1) # THIS SHOULD GIVE 1, BUT IT DOES NOT!!
        print square_average
        return 1 - average*average
 