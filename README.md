MPS on Probability
===

This program studies the stochastic process problems with an algorithm called the matrix product state (MPS) that showed its power first in many body quantum physics. Noticing the analogy between probability theory and quantum mechanics, we adapted the algorithm to deal with many variable stochastic process problems.

We can use transitional matrix methods to study stochastic processes and extract a lot of interesting properties. However, just like in the many body quantum system, the transitional matrix grows exponentially in a system of multiple degrees of freedom and becomes intractible with 30~50 variables. However, many models benefit from the property of locality and it is proved that, writing states as product of small matrices can effectively reduce the system complexity. We write transitional matrices as well as state vectors in terms of product of matrices and solve the time evolution problem efficiently.

A detailed description of the problems, models and the algorithm is in the documentation file Project_note.tex and a detailed user's manual will be available before the mid-January, 2015. If you want to play with the program before the beta version is released. 

Here is a quick start guide.
---
* proba.py is the example main caller and you can tune the parameters in the main function. It runs both exact and mps models. It calculates the joint probability, the mean value and the variance in both models.
* model.py defines pure model class for all the model classes. 
* exactsolver.py and mpssolver.py define two solvers: the former is the traditional transitional matrix method and the latter is the MPS solver. You can compare them for small sizes.

* measurement.py defines the base measurement classe Measurement for all the measurement classes. It can compute the joint probability, correlation function, variance, and mean value at any time point. For the last two calculations, user can define specific values for being angry (up) or calm (down), and by default, they are 1 and -1 respectively. Both time and site number counts from 0, and a python-way of specifying the number from the end of list is also accepted. Each measurement task must be specified in a particular way, though measurement time can be omitted, in which case it is taken to be the most recent time computed. The class uses the solver to initialize, and the measure functions take the same input. They have different requested tasks and each measurement function should check the task properties. It will raise an exception if it does not fit. 
  For example, ("correlation", 1, [1, 2, 3, 4, 5]) measures <S1*S2*S3*S4*S5> at t = 1, task_type_string = "Correlation", "Mean", "Variance", "Proba", more to be added some examples:("Correlation", 1, [1, 2, 3, 4, 5]) ("Mean", 1, [1]) ("Variance", 1, [1])("Proba", [(1, "up"), (3, "down")])It contains the following functions:
    1 addMeasureTask(self, task): The user can add different tasks by this function. The task is a tuple that shows up in the form (task_type_string, time_point, [points of interest]). Note in each task tuple, the integer following the task name will always be taken as the time step for all measurement. If time_point is not given, then it is assumed to be the last time step in calculation. It can also be negative which then follows Python's convention.
    2 measureProba(self, task): This function implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
    3 measureCorrelation(self, task, up=None, down=None): This implements the measurement of the n-point correlation function
    4 measureMean(self, task, up=None, down=None): It calculates the mean value of the model.
    5 measureVariance(self, task, up=None, down=None): Variance is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation.
    6 getTimeTask(self, task, task_temp): For a given task, find the time for the computation and the true tasks. If the second element in the task is an integer, then it is taken as the time. If time is not given in the task at the second position, then it is the last one. task_temp is the list that stores true tasks.



* exactmeasurement.py and mpsmeasurement.py are derived measurement classes that can compute mean, variance, correlation functions and joint probability for each model.
* proba_mps.py is the caller 

This project is still in progress and is subject to the GPL 2 licence. Please feel free to contact the authors if you have some questions or suggestions.

Authors:
---
* Bin Xu (binarybin)
* Peiqi Wang (pqwang1026)
* Liangsheng Zhang (phzls)
* Peter Gjeltema (PJ)
* Jun Xiong (xiongPU)
