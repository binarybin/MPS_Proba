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
        #check for correct task
        if task[0] != "Proba":
            raise Exception("Task is wrong type")
        #get tasks
        temptask = []
        time = self.getTimeTask(task, temptask)

        #negative sites case
        for i in range(len(temptask)):
            if temptask[i][0] < 0:
                temptask[i] = (self.solver.size + temptask[i][0], temptask[i][1])
        
        #compute probabilities
        proba = 0
        for i in range(len(self.basis)):
            match = True
            for j in range(len(temptask)):
                if self.basis[i][temptask[j][0]-1] != self.convert(temptask[j][1]):
                    match = False
            if match:
                proba += self.solver.results[time][i]
        return proba

    def measureCorrelation(self, task,up=None,down=None):
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

        for i in range(len(temptask)):
            if temptask[i] < 0:
                temptask[i] = self.solver.size + temptask[i]
        
        corr = 1
        task_index = 0
        for i in range(len(self.basis)):
            match = True
            for j in range(len(temptask)):
                if self.basis[i][temptask[j][0]-1] != self.convert(temptask[j][1]):
                    continue
                else:
                    corr *= self.solver.results[time][i]
        return corr
                
    def measureMean(self, task, up=None,down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Mean":
            raise Exception("Task is wrong type")

        if up is None:
            up = 1
        if down is None:
            down = -1

        temptask = []
        time = self.getTimeTask(task, temptask)
        newtask = ("", time, temptask)
        return self.measureCorrelation(newtask, up, down)
        
        
    def measureVariance(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        if task[0] != "Variance":
            raise Exception("Task is wrong type")
    
        if up is None:
            up = 1
        if down is None:
            down = -1

        temptask = []
        time = self.getTimeTask(task, temptask)
        newtask = ("", time, temptask)
        average = self.measureMean(newtask, up, down)
        up *= up
        down *= down
        square_average = self.measureMean(newtask, up, down)
        return square_average - average*average