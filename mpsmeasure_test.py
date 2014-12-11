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
            mean = 1-2*p

            task_var = ("Variance", [0])
            var = 1 - (1-2*p)**2

            task = [task_up, task_down, task_corr, task_mean, task_var]
            result = [prob_up, prob_down, corr, mean, var]

            for n in task:
                measure.addMeasureTask(n)

            measure.measure()

            for i in xrange(len(result)):
                self.assertAlmostEqual(measure.measure_result_list[i], result[i])

    def testTwosties(self):
        solver = TestSolver(2)
        output1 = []
        output2 = []

        w1 = np.ndarray(shape = (2, 2), dtype = float)
        w2 = np.ndarray(shape = (2, 2), dtype = float)

        for p1 in np.linspace(0,1,10):
            w1[0] = [p1/2, p1]
            w1[1] =  [(1-p1)/2, 1-p1]

            for p2 in np.linspace(0,1,10):
                w2[0] = [p2, p2/2]
                w2[1] =  [1-p2, (1-p2)/2]

                solver.mps_result = [[],[w1,w2]]

                measure = MpsMeasurement(solver, output1, output2)

                task_up_up = ("Proba", [(0, "up"), (1, "up")])
                prob_up_up = (1-p1) * (1-p2)

                task_down_up = ("Proba", [(0, "down"), (1, "up")])
                prob_down_up = p1 * (1-p2)

                task_corr = ("Correlation", [0,1])
                corr = (1-p1)*(1-p2) + p1*p2 - (1-p1)*p2 - (1-p2)*p1

                task_mean = ("Mean", [0])
                mean = 1 - 2*p1

                task_var = ("Variance", [1])
                var = 1 - (1-2*p2)**2

                task = [task_up_up, task_down_up, task_corr, task_mean, task_var]
                result = [prob_up_up, prob_down_up, corr, mean, var]

                for n in task:
                    measure.addMeasureTask(n)

                measure.measure()

                for i in xrange(len(result)):
                    self.assertAlmostEqual(measure.measure_result_list[i], result[i])



if __name__ == '__main__':
    unittest.main()
