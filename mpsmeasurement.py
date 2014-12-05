"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: mpsmeasurement.py
File purpose: the class for measuring several quantities using the exact solution
Responsible person: Liangsheng Zhang
"""

from measurement import Measurement
from exactsolver import ExactSolver
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

        up = 1
        down = 0

        # If time is not given in the task at the second position, then it is the last one
        if len(task) == 3:
            time = task[1]

            if time > len(self.solver.mps_result)-1 or -time > len(self.solver.mps_result):
                raise Exception("Time for probability calculation too large.")
        else:
            time = -1

        # Sort the task in ascending order of positions
        task_sort = sorted(task[1], key = lambda x:x[0])

        # Convert string to numbers
        for i in xrange(len(task_sort)):
            task_sort[i][1] = eval(task_sort[i][1])

        # Assume every site has same physical dimension d, then all possible physical states can
        # be expressed as a number written in d basis
        d = self.solver.mps_result[time][0].shape[0]
        permu_max = d**(self.solver.model.size - len(task[1]))

        # This number encodes all possible scenarios of physical indices
        permu = 0

        # The probability that would be returned
        prob = 0

        while permu < permu_max:
            task_index = 0
            prob_m = 1 # A temperory val storing the prob of a given configuration
                       # indicated by permu

            permu_temp = permu

            for i in range(self.solver.model.size):
                mps_temp = self.solver.mps_result[time][i]

                # This is the site mentioned in the task
                if task_index < len(task_sort) and i == task_sort[task_index]:
                    spin = task_sort[task_index][1]
                    prob_m = np.dot(prob_m, mps_temp[spin])
                    task_index += 1

                else:
                    r = permu_temp % d # The physical index given by permu at this site
                    permu_temp /= d
                    prob_m = np.dot(prob_m, mps_temp[r])

            permu += 1
            prob += prob_m

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

        time = task[1]
        sites = task[2]

        # If time is not given in the task at the second position, then it is the last one
        if len(task) == 3:
            time = task[1]

            if time > len(self.solver.mps_result)-1 or -time > len(self.solver.mps_result):
                raise Exception("Time for correlation calculation too large.")
        else:
            time = -1

        # Assume every site has same physical dimension d, then all possible physical states can
        # be expressed as a number written in d basis
        d = self.solver.mps_result[time][0].shape[0]
        permu_max = d**(self.solver.model.size)

        # This number encodes all possible scenarios of physical indices
        permu = 0

        # The probability that would be returned
        corr = 0

        while permu < permu_max:
            task_index = 0
            prob_m = 1 # A temperory val storing the prob of a given configuration
                       # indicated by permu

            permu_temp = permu
            val = 1

            for i in range(self.solver.model.size):
                mps_temp = self.solver.mps_result[time][i]

                r = permu_temp % d # The physical index given by permu at this site
                permu_temp /= d
                prob_m = np.dot(prob_m, mps_temp[r])

                # This is the site mentioned in the task
                if task_index < len(sites) and i == sites[task_index]:
                    val *= r*up + (1-r)*down
                    task_index += 1

            permu += 1
            corr += prob_m * val

        return corr

    def measureMean(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        # The value of up and down spins
        if up is None:
            up = 1
        if down is None:
            down = -1

        # If time is not given in the task at the second position, then it is the last one
        if len(task) == 3:
            time = task[1]

            if time > len(self.solver.mps_result)-1 or -time > len(self.solver.mps_result):
                raise Exception("Time for mean calculation too large.")
        else:
            time = -1

        task_new = ("", time, range(self.solver.size))

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

        # If time is not given in the task at the second position, then it is the last one
        if len(task) == 3:
            time = task[1]

            if time > len(self.solver.mps_result)-1 or -time > len(self.solver.mps_result):
                raise Exception("Time for mean calculation too large.")
        else:
            time = -1

        task_new = ("", time)
        ave = self.measureMean(task_new, up, down)

        up *= up
        down *= down

        square_ave = self.measureMean(task_new, up, down)

        return square_ave - ave*ave
