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
from angryboys import AngryBoys
from radiatingboys import RadiatingBoys
from exponentialboys import ExponentialBoys
from exactsolver import ExactSolver
from mpssolver import MpsSolver
from exactmeasurement import ExactMeasurement
from mpsmeasurement import MpsMeasurement
import pylab

if __name__=="__main__":
    """
    The main call of this program.
    """
    total_time = 20
    main_output = sys.stdout # the output channel for important results that will be shown in the final run
    aux_output = sys.stdout # the output channel for auxiliary information that may be interesting for dubugging but not in the final run

    angry_boys = AngryBoys(size = 10, remain_proba = 0.1, init_state = "all down", output1 = main_output, output2 = aux_output)
#    angry_boys = RadiatingBoys(size = 10, remain_proba = 0.1, nearest_neighbour_proba = 0.4, second_neighbour_proba = 0.5, init_state = "all down", output1 = main_output, output2 = aux_output)
#    angry_boys = ExponentialBoys(size = 10, J = 0.5, K = 1, init_state = "all down", output1 = main_output, output2 = aux_output)
    exact_solver = ExactSolver(model = angry_boys, output1 = main_output, output2 = aux_output)
    mps_solver = MpsSolver(model = angry_boys, bound_dimension = 5, output1 = main_output, output2 = aux_output)
    exact_solver.evolve(total_time)
    mps_solver.evolve(total_time)

    exact_measurement = ExactMeasurement(exact_solver, main_output, aux_output)
    mps_measurement   = MpsMeasurement(mps_solver, main_output, aux_output)

    for i in range(total_time):
        exact_measurement.addMeasureTask(("Proba", i, [(5, "up")]))
        mps_measurement.addMeasureTask(("Proba", i, [(5, "up")]))

    exact_measurement.measure()
    mps_measurement.measure()

    ers = exact_measurement.measure_result_list
    mrs = mps_measurement.measure_result_list

    pylab.plot(range(len(ers)), ers, label = "exact")
    pylab.plot(range(len(mrs)), mrs, label = "mps")
    pylab.legend(loc="lower right")
    pylab.show()