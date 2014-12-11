from mpsmeasurement import MpsMeasurement
import unittest
import numpy as np

class TestSolver(object):
    def __init__(self, L):
        self.L = L

class MpsMeasurementTest(unittest.TestCase):
    def testOnestie(self):
        solver = TestSolver(1)
        output1 = []
        output2 = []

        w = np.ndarray(shape = (2, 1), dtype = float)

        for p in np.linspace(0,1,10):
            w[0] = p
            w[1] = 1-p

            solver.mps_result = [[]]
            solver.mps_result[0].append(w)

            measure = MpsMeasurement(solver, output1, output2)

            task_up = ("Proba", (0, "up"))
            prob_up = 1-p

            task_down = ("Proba", (0, "down"))
            prob_down = p

            task_corr = ("Correlation", [0])
            corr = 1-2*p

            task_mean = ("Mean", [0])
            mean = measure.measureCorrelation(task_mean)

            task_var = ("Variance", [0])
            var = 1 - (1-2*p)**2

            task = [task_up, task_down, task_corr, task_mean, task_var]
            result = [prob_up, prob_down, corr, mean, var]

            for n in task:
                measure.addMeasureTask(n)

            measure.measure()

            for i in xrange(len(result)):
                self.assertAlmostEqual(measure.measure_result_list[i], result[i])



if __name__ == '__main__':
    unittest.main()
