"""
Program name: MPS-Proba
Test purpose: Plot joint proba v.s. time at different bound dimensions
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
import csv


if __name__=="__main__":
    """
    The main call of this program.
    """
    total_time = 100

#    angry_boys = AngryBoys(size = 200, remain_proba = 0.1, init_state = "all down")
#    angry_boys = RadiatingBoys(size = 200, remain_proba = 0.1, nearest_neighbour_proba = 0.4, second_neighbour_proba = 0.5, init_state = "all down")
    angry_boys = ExponentialBoys(size = 40, J = 0.5, K = 0.5, init_state = "all down")
#    angry_boys = ProjectionBoys(size = 100, p0 = 1.0, p1 = 2.0, q1 = 2.0, q2 = 2.0, init_state = "all down")

    mps_solver = MpsSolver(model = angry_boys, bound_dimension = 10)
    mps_solver.evolve(total_time)
    mps_measurement   = MpsMeasurement(mps_solver)
    mps_solver2 = MpsSolver(model = angry_boys, bound_dimension = 20)
    mps_solver2.evolve(total_time)
    mps_measurement2   = MpsMeasurement(mps_solver2)
#    mps_solver3 = MpsSolver(model = angry_boys, bound_dimension = 5)
#    mps_solver3.evolve(total_time)
#    mps_measurement3   = MpsMeasurement(mps_solver3)
#    mps_solver4 = MpsSolver(model = angry_boys, bound_dimension = 2)
#    mps_solver4.evolve(total_time)
#    mps_measurement4   = MpsMeasurement(mps_solver4)


    for i in range(total_time):
        mps_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
        mps_measurement2.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
#        mps_measurement3.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
#        mps_measurement4.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
    mps_measurement.measure()
    mpb = deepcopy(mps_measurement.measure_result_list)
    mps_measurement.clearMeasurement()

    mps_measurement2.measure()
    mpb2 = deepcopy(mps_measurement2.measure_result_list)
    mps_measurement2.clearMeasurement()

#    mps_measurement3.measure()
#    mpb3 = deepcopy(mps_measurement3.measure_result_list)
#    mps_measurement3.clearMeasurement()
#
#    mps_measurement4.measure()
#    mpb4 = deepcopy(mps_measurement4.measure_result_list)
#    mps_measurement4.clearMeasurement()

#    for i in range(total_time):
#        mps_measurement.addMeasureTask(("Mean", i, [5]))
#    mps_measurement.measure()
#    mmn = deepcopy(mps_measurement.measure_result_list)
#    mps_measurement.clearMeasurement()

#    for i in range(total_time):
#        mps_measurement.addMeasureTask(("Variance", i, [5]))
#    mps_measurement.measure()
#    mva = deepcopy(mps_measurement.measure_result_list)
#    mps_measurement.clearMeasurement()
#    mpb2a=[j for j in mpb]

    with open("Data/Exponential_MPS_t100_s40_bd10to20b.csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerows([mpb,mpb2])
    pylab.plot(range(len(mpb)), mpb, "-o", label = "bound dimension = 10")
    pylab.plot(range(len(mpb2)), mpb2, "-*", label = "bound dimension = 20")
#    pylab.plot(range(len(mpb)), mpb3, "-s", label = "bound dimension = 5")
#    pylab.plot(range(len(mpb2)), mpb4, "-^", label = "bound dimension = 2")
#    pylab.plot(range(len(mmn)), mmn, "-o", label = "mps, mean")
#    pylab.plot(range(len(mva)), mva, "-o", label = "mps, variance")

    pylab.legend(loc="lower right")
    pylab.xlabel("time")
    pylab.ylabel("joint probability")
    pylab.ylim([0,0.35])
    pylab.savefig('Figs/Exponential_MPS_t100_s40_bd10to20b.eps')
    pylab.show()

