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
from projectionboys import ProjectionBoys
import pylab

if __name__=="__main__":
    """
    The main call of this program.
    """
    total_time = 40

#    angry_boys = AngryBoys(size = 10, remain_proba = 0.1, init_state = "all down")
    angry_boys = RadiatingBoys(size = 10, remain_proba = 0.1, nearest_neighbour_proba = 0.4, second_neighbour_proba = 0.5, init_state = "all down")
#    angry_boys = ExponentialBoys(size = 10, J = 0.5, K = 0.5, init_state = "all down")
#    angry_boys = ProjectionBoys(size = 100, p0 = 1.0, p1 = 2.0, q1 = 2.0, q2 = 2.0, init_state = "all down")

    """
    mrs = []
    evolve_time = []

    import time
    run_list = range(3,8)
    for i in run_list:
        start_time = time.time()
        print "bound dimension: ", i
        print "Solver Initialized"
        mps_solver = MpsSolver(model = angry_boys, bound_dimension = i)
        print "Evolve"
        mps_solver.evolve(total_time)
        print "Measure"
        mps_measurement = MpsMeasurement(mps_solver)
        mps_measurement.addMeasureTask(("Proba", [(5, "up"), (-1, "down")]))
        mps_measurement.measure()
        mrs.append(mps_measurement.measure_result_list[0])
        evolve_time.append(time.time()-start_time)
        print "Execution time: ", evolve_time[-1]

    pylab.figure(1)
    pylab.plot(run_list, mrs)

    pylab.figure(2)
    pylab.plot(run_list,evolve_time)
    pylab.show()

    """
    """
    for d in [12, 13, 14, 15, 16]:
        mps_solver = MpsSolver(model = angry_boys, bound_dimension = d)
        mps_solver.evolve(total_time)
        mps_measurement   = MpsMeasurement(mps_solver)
        for i in range(total_time):
            mps_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
        mps_measurement.measure()
        mrs = mps_measurement.measure_result_list
        pylab.plot(range(len(mrs)), mrs, label = "mps d = "+str(d))
    
    
    """
    exact_solver = ExactSolver(model = angry_boys)
    mps_solver = MpsSolver(model = angry_boys, bound_dimension = 10)

    exact_solver.evolve(total_time)
    mps_solver.evolve(total_time)

    exact_measurement = ExactMeasurement(exact_solver)
    mps_measurement   = MpsMeasurement(mps_solver)


    for i in range(20):
    #    exact_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1,"down")]))
    #    mps_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
        exact_measurement.addMeasureTask(("Variance", i, [3]))
        mps_measurement.addMeasureTask(("Variance", i, [3]))

    exact_measurement.measure()
    mps_measurement.measure()

    ers = exact_measurement.measure_result_list
    mrs = mps_measurement.measure_result_list

    pylab.plot(range(len(ers)), ers, "-", label = "exact")
    pylab.plot(range(len(mrs)), mrs, "o", label = "mps")
    
    pylab.legend(loc="lower right")
    pylab.show()
    

