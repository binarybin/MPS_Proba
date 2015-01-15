"""
    Program name: MPS-Proba
    Program purpose: Used to calculate the error between Exact and MPS. The Alpha version of the APC 524 project.
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
import numpy as np

if __name__=="__main__":
    """
        The main call of this program.
        """
    total_time = 100
    
    angry_boys = AngryBoys(size = 10, remain_proba = 0.1, init_state = "all down")
#    angry_boys = RadiatingBoys(size = 10, remain_proba = 0.1, nearest_neighbour_proba = 0.4,second_neighbour_proba = 0.5, init_state = "all down")
#    angry_boys = ExponentialBoys(size = 10, J = 0.5, K = 0.5, init_state = "all down")
#    angry_boys = ProjectionBoys(size = 10, p0 = 1.0, p1 = 2.0, q1 = 2.0, q2 = 2.0, init_state = "all down")

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
    prob_err_bd=[]
    mean_err_bd=[]
    var_err_bd=[]
    for bd in range(1,11):
        prob_err_sum=0
        mean_err_sum=0
        var_err_sum=0
        exact_solver = ExactSolver(model = angry_boys)
        mps_solver = MpsSolver(model = angry_boys, bound_dimension = bd+1)
        
        exact_solver.evolve(total_time)
        mps_solver.evolve(total_time)
        
        exact_measurement = ExactMeasurement(exact_solver)
        mps_measurement   = MpsMeasurement(mps_solver)
        
        
        for i in range(total_time):
            exact_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1,"down")]))
            mps_measurement.addMeasureTask(("Proba", i, [(5, "up"), (-1, "down")]))
        exact_measurement.measure()
        mps_measurement.measure()

        epb = deepcopy(exact_measurement.measure_result_list)
        mpb = deepcopy(mps_measurement.measure_result_list)
        exact_measurement.measure_result_list = []
        mps_measurement.measure_result_list = []
        exact_measurement.measurement_list = []
        mps_measurement.measurement_list = []
        
        for i in range(total_time):
            exact_measurement.addMeasureTask(("Mean", i, [5]))
            mps_measurement.addMeasureTask(("Mean", i, [5]))
        exact_measurement.measure()
        mps_measurement.measure()

        emn = deepcopy(exact_measurement.measure_result_list)
        mmn = deepcopy(mps_measurement.measure_result_list)
        exact_measurement.measure_result_list = []
        mps_measurement.measure_result_list = []
        exact_measurement.measurement_list = []
        mps_measurement.measurement_list = []
        
        for i in range(total_time):
            exact_measurement.addMeasureTask(("Variance", i, [5]))
            mps_measurement.addMeasureTask(("Variance", i, [5]))
        exact_measurement.measure()
        mps_measurement.measure()

        eva = deepcopy(exact_measurement.measure_result_list)
        mva = deepcopy(mps_measurement.measure_result_list)
        exact_measurement.measure_result_list = []
        mps_measurement.measure_result_list = []
        exact_measurement.measurement_list = []
        mps_measurement.measurement_list = []
        
        mpb_arr=np.asarray(mpb)
        epb_arr=np.asarray(epb)
        mmn_arr=np.asarray(mmn)
        emn_arr=np.asarray(emn)
        mva_arr=np.asarray(mva)
        eva_arr=np.asarray(eva)
        prob_err_sum=sum((mpb_arr - epb_arr)**2)
        mean_err_sum=sum((mmn_arr - emn_arr)**2)
        var_err_sum=sum((mva_arr - eva_arr)**2)
        prob_err_bd.append(prob_err_sum)
        mean_err_bd.append(mean_err_sum)
        var_err_bd.append(mean_err_sum)

#        prob_err_bd=np.append(prob_err_bd, prob_err_sum)
#        mean_err_bd=np.append(mean_err_bd, mean_err_sum)
#        var_err_bd=np.append(var_err_bd, mean_err_sum)


    pylab.plot(range(1,(len(prob_err_bd)+1)), prob_err_bd, "-o", label = "joint probability square error sum")
#    pylab.plot(range(1,(len(mean_err_bd)+1)), mean_err_bd, "o", label = "mean square error sum")
#    pylab.plot(range(1,(len(var_err_bd)+1)), var_err_bd, "o", label = "variance square error sum")

    
    pylab.legend(loc="lower left")
    pylab.xlabel("Bound Dimension")
    pylab.ylabel("Square Error Sum")
    pylab.yscale('log')
#    pylab.yticks(np.asarray([1E-6,1E-5,1E-4,1E-3,1E-2]))
    pylab.savefig('Figs/Angry_Error_t100_s10_bd10_log.eps')
    #pylab.savefig('Figs/1.eps')
    pylab.show()