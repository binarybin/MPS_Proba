"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: measurement.py
File purpose: the abstract measurement class
Responsible person: Bin Xu
"""
from exactsolver import ExactSolver
from mpssolver import MpsSolver

class Measurement(object):
    """
    The abstract Measurement base class
    """

    def __init__(self, solver, output1, output2):
        self.solver = solver
        self.output1 = output1
        self.output2 = output2
        self.measurement_list = []
        self.measure_result_list = []

    def addMeasureTask(self, task):
        self.measurement_list.append(task)

    def measure(self):
        for task in self.measurement_list:
            self.measure_result_list.append(eval("self.measure"+task[0]+"(task)"))


    # All measure functions take the same input. They have different requested tasks and each measurement function should check the task properties. Raise an exception if
    # it does not fit.

    # task is a tuple that shows up in the form (task_type_string, time_point, [points of interest]), if time_point is not given, then it is assumed to be the last time step in calculation. It can also be negative which then follows Python's convention.
    # for example, ("correlation", 1, [1, 2, 3, 4, 5]) measures <S1*S2*S3*S4*S5> at t = 1
    # task_type_string = "Correlation", "Mean", "Variance", "Proba", more to be added
    # some examples:
    # ("Correlation", 1, [1, 2, 3, 4, 5])
    # ("Mean", 1, [1])
    # ("Variance", 1, [1])
    # ("Proba", [(1, "up"), (3, "down")])
    def measureProba(self, task):
        """
        This implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
        """
        raise NotImplementedError("Unimplemented abstract method")

    def measureCorrelation(self, task, up=None, down=None):
        """
        This implements the measurement of the n-point correlation function
        """
        raise NotImplementedError("Unimplemented abstract method")

    def measureMean(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        raise NotImplementedError("Unimplemented abstract method")

    def measureVariance(self, task, up=None, down=None):
        """
        Mean is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation
        """
        raise NotImplementedError("Unimplemented abstract method")