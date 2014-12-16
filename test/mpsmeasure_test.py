from mpsmeasurement import MpsMeasurement
import unittest
import numpy as np
from math import sqrt

class TestSolver(object):
    def __init__(self, L):
        self.L = L

class MpsMeasurementTest(unittest.TestCase):
    def testOnestie(self):
        solver = TestSolver(1)
        output1 = []
        output2 = []

        w = np.ndarray(shape = (2, 1, 1), dtype = float)

        for p in np.linspace(0,1,10):
            w[0] = [p]
            w[1] = [1-p]

            solver.results = [[]]
            solver.results[0].append(w)

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

        w1 = np.ndarray(shape = (2, 1, 2), dtype = float)
        w2 = np.ndarray(shape = (2, 2, 1), dtype = float)

        for p1 in np.linspace(0,1,10):
            w1[0] = [[p1/2, p1]]
            w1[1] = [[(1-p1)/2, 1-p1]]

            for p2 in np.linspace(0,1,10):
                w2[0] = [[p2], [p2/2]]
                w2[1] = [[1-p2], [(1-p2)/2]]

                solver.results = [[w1,w2]]

                measure = MpsMeasurement(solver, output1, output2)

                task_up_up = ("Proba", (0, "up"), (1, "up"))
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

    def testThreesites(self):
        solver = TestSolver(3)
        output1 = []
        output2 = []

        w1 = np.ndarray(shape = (2, 1, 2), dtype = float)
        w2 = np.ndarray(shape = (2, 2, 4), dtype = float)
        w3 = np.ndarray(shape = (2, 4, 1), dtype = float)

        w1[0] = [[0,1]]
        w1[1] = [[1,0]]

        for p in np.linspace(0.1,0.5,10):
            for q in np.linspace(0.1,0.5,10):
                r = (1 - p -q)/2

                w2[0] = ([[q/sqrt(p**2 + q**2), 0, 0, p/sqrt(p**2 + q**2)],
                         [0, 1.0/sqrt(2), 1.0/sqrt(2), 0]])
                w2[1] = ([[p/sqrt(p**2 + q**2), 0, 0, -p/sqrt(p**2+q**2)],
                         [0, 1.0/sqrt(2), -1.0/sqrt(2), 0]])
                w3[0] = [[0], [sqrt(2)*r], [0], [0]]
                w3[1] = [[sqrt(p**2 + q**2)], [0], [0], [0]]

                solver.results = [[w1, w2, w3]]

                measure = MpsMeasurement(solver, output1, output2)

                task_up = ("Proba", (0, "up"))
                prob_up = p + q

                task_r_up_down = ("Proba", (1, "up"), (-1, "down"))
                prob_r_up_down = r

                task_corr1 = ("Correlation",[0, 2])
                corr1 = 1

                task_corr2 = ("Correlation", [1,-1])
                corr2 = p - q

                task_mean1 = ("Mean", [1])
                mean1 = p - q

                task_mean2 = ("Mean", [0])
                mean2 = p + q - 2 * r

                task_var = ("Variance", [-1])
                var = 1 - (p+q-2*r)**2

                task = ([task_up, task_r_up_down, task_corr1, task_corr2, task_mean1,
                         task_mean2, task_var])
                result = [prob_up, prob_r_up_down, corr1, corr2, mean1, mean2, var]

                for n in task:
                    measure.addMeasureTask(n)

                measure.measure()

                for i in xrange(len(result)):
                    self.assertAlmostEqual(measure.measure_result_list[i], result[i])


if __name__ == '__main__':
    unittest.main()
