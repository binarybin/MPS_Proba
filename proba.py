"""
Program name: MPS-Proba
Program purpose: The Alpha version of the APC 524 project.
Due date: 12/12/2014
Author list:
    Bin Xu
    Liangsheng Zhang
    Peiqi Wang
    Peter Gjeltema
    Jun Xiong
File list(author initial):
    proba.py -- the main routine that wraps everything (BX)
    exact.py -- the exact solution class (PG)
    compression.py -- the MPS compression class (LZ & JX)
    measurement.py -- the measurement of some interesting quantities (PW)

Code style conventions:
    variable_name
    functionName
    ClassName
    CONSTANT_NAME
    These standard should be followed strictly, for example, a class is called MpsSolver (not MPSSolver) and an instance of it is called mps_solver (not MPS_solver)
    
    Please define the interface functions in a way that they take all variables with explicit name
"""

import sys

if __name__=="main":
    """
    The main call of this program.
    """
    
    main_output = sys.stdout # the output channel for important results that will be shown in the final run
    aux_output = sys.stdout  # the output channel for auxiliary information that may be interesting for dubugging but not in the final run
    
    angry_boy = ProbaModel(size = 10, change_rate = 0.1, init_state = "all down", output1 = main_output, output2 = aux_output)
    exact_solver = ExactSolver(model = angry_boy, output1 = main_output, output2 = aux_output)
    mps_solver = MpsSolver(model = angry_boy, bound_dimension = 20, output1 = main_output, output2 = aux_output)
    exact_solver.evolve(step = 20)
    mps_solver.evolve(step = 20)
    
    exact_measurement = ExactMeasurement(exact_solver, main_output, aux_output)
    mps_measurement   = MpsMeasurement(mps_solver, main_output, aux_output)


    
