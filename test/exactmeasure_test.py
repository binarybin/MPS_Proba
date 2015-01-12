import sys
sys.path.append('../src/')

from exactmeasurement import ExactMeasurement
import unittest
import numpy as np

class TestModel(object):
    def __init__(self, L):
        self.size = L

class TestSolver(object):
    def __init__(self, L):
        self.model = TestModel(L)

class ExactMeasurementTest(unittest.TestCase):
    def testOnestie(self):
        solver = TestSolver(1)

        w = np.ndarray(shape = (2, 1, 1), dtype = float)

        for p in np.linspace(0,1,10):
            w = np.matrix([[p],[1-p]])

            solver.results = [w]

            measure = ExactMeasurement(solver)

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

        for p1 in np.linspace(0,1,10):
            for p2 in np.linspace(0,1,10):
                w = np.matrix([[p1*p2],[p1*(1-p2)],[(1-p1)*p2],[(1-p1)*(1-p2)]])

                solver.results = [w]

                measure = ExactMeasurement(solver)

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

        for p in np.linspace(0.1,0.5,10):
            for q in np.linspace(0.1,0.5,10):
                r = (1 - p -q)/2

                w = np.matrix([[r],[0],[r],[0],[0],[q],[0],[p]])

                solver.results = [w]

                measure = ExactMeasurement(solver)

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