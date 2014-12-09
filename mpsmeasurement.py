"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpsmeasurement.py
File purpose: the class for measuring several quantities using the exact solution
Responsible person: Liangsheng Zhang
"""

from measurement import Measurement
from exactsolver import ExactSolver
from mpssolver import MpsSolver
import numpy as np

class MpsMeasurement(Measurement):
    """
    The Measurement class for the mps method
    """
    def measureProba(self, task):
        """
        This implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
        """
        # Assume for a mps, physical index 0 represents down and physical index 1 represents up
        # Assume solver.mps_result has elements time ordered and start from time 0

        up = 1 # physical index for spin up
        down = 0 # physical index for spin down

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        # Sort the task in ascending order of positions
        task_sort = sorted(task_temp, key = lambda x:x[0])

        # Convert string to numbers
        for i in xrange(len(task_sort)):
            task_sort[i][1] = eval(task_sort[i][1])

        # Use distributive law to compute probability
        task_index = 0 # Position in task_sort          
        prob = 1 # The probability that would be returned

        for i in xrange(self.solver.L):
            mps_temp = self.solver.mps_result[time][i]
            if task_index < len(task_sort) and i == task_sort[task_index]:
                spin = task_sort[task_index][1]
                prob_mps = mps_temp[spin]
                task_index += 1
            else:
                prob_mps = mps_temp[0] + mps_temp[1]
            prob = np.dot(prob, prob_mps)

        return prob

    def measureCorrelation(self, task, up=None, down=None):
        """
        This implements the measurement of the n-point correlation function. User can specify
        the values corresponding to up and down spins. If they are not given, then by default
        up is 1 and down is -1
        """
        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        sites = task_temp.sort()

        # Use distributive law to compute probability
        task_index = 0 # Position in sites          
        corr = 1 # The correlation that would be returned

        for i in xrange(self.solver.L):
            mps_temp = self.solver.mps_result[time][i]
            if task_index < len(sites) and i == sites[task_index]:
                prob_mps = down * mps_temp[0] + up * mps_temp[1]
                task_index += 1
            else:
                prob_mps = mps_temp[0] + mps_temp[1]
            corr = np.dot(corr, prob_mps)

        return corr

    def measureMean(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type
        and call measureCorrelation
        """
        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        task_new = ("", time, task_temp)

        return self.measureCorrelation(task_new, up, down)

    def measureVariance(self, task, up=None, down=None):
        """
        measureMean is used when evaluating variance.
        """

        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        task_new = ("", time, task_temp)

        ave = self.measureMean(task_new, up, down)

        up *= up
        down *= down

        square_ave = self.measureMean(task_new, up, down)

        return square_ave - ave*ave

    def getTimeTask(self, task, task_temp):
        """ 
        For a given task, find the time for the computation and the true tasks.
        If time is not given in the task at the second position, then it is the last one.
        task_temp is the list that stores true tasks.
        """

        if len(task) == 3:
            time = task[1]
            task_pos = 2 # The positive in task where sites are given

            if time > len(self.solver.mps_result)-1 or -time > len(self.solver.mps_result):
                raise Exception("Time for calculation is too large.")
        else:
            time = -1
            task_pos = 1

        if type(task[task_pos]) is list:
            task_temp = task[task_pos]
        else:
            while task_pos < len(task):
                task_temp.append(task[task_pos])
                task_pos += 1

        return time