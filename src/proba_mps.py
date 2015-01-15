"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
File name: proba.py
File purpose: the main routine
Responsible person: Bin Xu
Due date: 12/12/2014
Author list:
    Bin Xu
    Liangsheng Zhang
    Peiqi Wang
    Peter Gjeltema
    Jun Xiong
File list(author initial):
    proba.py -- the main routine that wraps everything
    exact.py -- the exact solution class
    compression.py -- the MPS compression class
    measurement.py -- the measurement of some interesting quantities

Code style conventions:
    variable_name
    functionName
    ClassName
    CONSTANT_NAME
    These standard should be followed strictly, for example, a class is called MpsSolver (not MPSSolver) and an instance of it is called mps_solver (not MPS_solver).
    Another example, a variable is called hamiltonian not Hamiltonian.

    Please define the interface functions in a way that they take all variables with explicit name
"""

import sys
from copy import deepcopy
from angryboys import AngryBoys
from radiatingboys import RadiatingBoys
from exponentialboys import ExponentialBoys
from exactsolver import ExactSolver
from mpssolver import MpsSolver
from exactmeasurement import ExactMeasurement
from mpsmeasurement import MpsMeasurement
from projectionboys import ProjectionBoys
import pylab

if __name__=="__main__":
    """
    The main call of this program.
    """
    total_time = 100

    angry_boys = AngryBoys(size = 200, remain_proba = 0.1, init_state = "all down")
#    angry_boys = RadiatingBoys(size = 200, remain_proba = 0.1, nearest_neighbour_proba = 0.4, second_neighbour_proba = 0.5, init_state = "all down")
#    angry_boys = ExponentialBoys(size = 200, J = 0.5, K = 0.5, init_state = "all down")
#    angry_boys = ProjectionBoys(size = 200, p0 = 1.0, p1 = 2.0, q1 = 2.0, q2 = 2.0, init_state = "all down")

    mps_solver = MpsSolver(model = angry_boys, bound_dimension = 10)
    mps_solver.evolve(total_time)
    mps_measurement   = MpsMeasurement(mps_solver)


    for i in range(total_time):
        mps_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
    mps_measurement.measure()
    mpb = deepcopy(mps_measurement.measure_result_list)
    mps_measurement.clearMeasurement()
    """
    for i in range(total_time):
        mps_measurement.addMeasureTask(("Mean", i, [5]))
    mps_measurement.measure()
    mmn = deepcopy(mps_measurement.measure_result_list)
    mps_measurement.clearMeasurement()
    
    for i in range(total_time):
        mps_measurement.addMeasureTask(("Variance", i, [5]))
    mps_measurement.measure()
    mva = deepcopy(mps_measurement.measure_result_list)
    mps_measurement.clearMeasurement()
    """
    pylab.plot(range(len(mpb)), mpb, "-o", label = "mps, joint probability")
#    pylab.plot(range(len(mmn)), mmn, "-o", label = "mps, mean")
 #   pylab.plot(range(len(mva)), mva, "-o", label = "mps, variance")
    
    pylab.legend(loc="lower right")
    pylab.xlabel("time")
    pylab.ylabel("Measurements")
    pylab.show()
