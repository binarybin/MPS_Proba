"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpsmeasurement.py
File purpose: the class for measuring several quantities using the exact solution
Responsible person: Liangsheng Zhang
"""

from measurement import Measurement
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
        # Assume solver.results has elements time ordered and start from time 0
        # Assume the sites start from 0
        # Assume the first mps is of shape (2,1,n)
        # Assume the last mps if of shape (2,n,1)

        if task[0] != "Proba":
            raise Exception("Proba function is called incorrectly in MPS measurement")

        up = 1 # physical index for spin up
        down = 0 # physical index for spin down

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        # Allow negative site number following python convention
        for index in range(len(task_temp)):
            if task_temp[index][0] < 0:
                task_temp[index] = (self.solver.L + task_temp[index][0], task_temp[index][1])


        # Sort the task in ascending order of positions
        task_sort = sorted(task_temp, key = lambda x:x[0])

        # Use distributive law to compute probability
        task_index = 0 # Position in task_sort
        prob = 1 # The probability that would be returned

        for i in xrange(self.solver.L):
            mps_temp = self.solver.results[time][i]
            if task_index < len(task_sort) and i == task_sort[task_index][0]:
                spin = eval(task_sort[task_index][1])
                prob_mps = mps_temp[spin]
                task_index += 1
            else:
                prob_mps = mps_temp[0] + mps_temp[1]
            prob = np.dot(prob, prob_mps)

        return float(prob)

    def measureCorrelation(self, task, up=None, down=None):
        """
        This implements the measurement of the n-point correlation function. User can specify
        the values corresponding to up and down spins. If they are not given, then by default
        up is 1 and down is -1
        """

        if task[0] != "Correlation":
            raise Exception("Correlation function is called incorrectly in MPS measurement")

        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        for index in range(len(task_temp)):
            if task_temp[index]<0:
                task_temp[index] = self.solver.L + task_temp[index]

        task_temp.sort()

        # Use distributive law to compute probability
        task_index = 0 # Position in sites
        corr = 1 # The correlation that would be returned

        for i in xrange(self.solver.L):
            mps_temp = self.solver.results[time][i]
            if task_index < len(task_temp) and i == task_temp[task_index]:
                prob_mps = down * mps_temp[0] + up * mps_temp[1]
                task_index += 1
            else:
                prob_mps = mps_temp[0] + mps_temp[1]
            corr = np.dot(corr, prob_mps)

        return float(corr)

    def measureMean(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type
        and call measureCorrelation
        """

        if task[0] != "Mean":
            raise Exception("Mean function is called incorrectly in MPS measurement")

        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        task_new = ("Correlation", time, task_temp)

        return self.measureCorrelation(task_new, up, down)

    def measureVariance(self, task, up=None, down=None):
        """
        measureMean is used when evaluating variance.
        """

        if task[0] != "Variance":
            raise Exception("Variance function is called incorrectly in MPS measurement")

        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # Get the time and true tasks
        task_temp = []
        time = self.getTimeTask(task, task_temp)

        task_new = ("Mean", time, task_temp)

        ave = self.measureMean(task_new, up, down)

        up *= up
        down *= down

        square_ave = self.measureMean(task_new, up, down)

        return square_ave - ave*ave
