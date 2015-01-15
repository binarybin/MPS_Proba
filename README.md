MPS on Probability
===

This program studies the stochastic process problems with an algorithm called the matrix product state (MPS) that showed its power first in many body quantum physics. Noticing the analogy between probability theory and quantum mechanics, we adapted the algorithm to deal with many variable stochastic process problems.

We can use transitional matrix methods to study stochastic processes and extract a lot of interesting properties. However, just like in the many body quantum system, the transitional matrix grows exponentially in a system of multiple degrees of freedom and becomes intractible with 30~50 variables. However, many models benefit from the property of locality and it is proved that, writing states as product of small matrices can effectively reduce the system complexity. We write transitional matrices as well as state vectors in terms of product of matrices and solve the time evolution problem efficiently.

A detailed description of the problems, models and the algorithm is in the documentation file Project_note.pdf and here is a quick start guide.

Structure:
---
* proba.py is an example main caller and you can tune the parameters in the main function which runs both exact and mps models. It calculates the joint probability, the mean value and the variance in different models at any time points. And it can compare the to exact solution for short chains. Please be aware that the exact method can only calculate a short chain of a typical maximum size of 14, and the program raises an exception if the chain is too long to solve exactly.
* proba_mps.py is an example that only uses the MPS method. Its function is similar to proba.py but do not perform an exact solution. You can feel free to try a system of size 200. It has four models and the user can choose which model to use. For example, angry_boys = AngryBoys(size = 200, remain_proba = 0.1, init_state = "all down") defines a chain using the Angryboys model with size 200.
* model.py is the base model class for all the model classes. 
* solver.py is the base solver class. 
* exactsolver.py defines the exact traditional matrix method.
* mpssolver.py defines the MPSsolver for the matrix product states method. It has the following variables that the user can tune according to his needs. 
L: length of the chain. 
bound_dimension: dimension of the auxilary space for the compressed MPS
n: dimension of the physical space
mps: current state expressed in mps. It is the mps obtained after apply the mpo
mpsc: current compressed mps. It is obtained by calling compression algorithm
mpo: operator. By convention we always apply mpo on mpsc
partial_overlap_lr, partial_overlap_rl: the partial overlap between mps and mpsc. this is only useful in the algorithm of compression by variation
results: list of mps keeping the history
t: current time
epsil: error threshold for CompressionVariational.
cpr_err: error of the compression(L2 distance between compressed state and true state).
The included functions are:
compressionVariational(self,direction=0,sweep_nbr=0,norm=0): compressesthe MPS stored in self.mps and writes the results in self.mpsc. The User can set the dimension of compressed MPS via the variable bound_dimension. The tolerance of convergence for sweeps can be accessed by self.epsil. Options available for self.compressionVariational include:
A. direction: the direction of first sweep. 0 for left to right and 1 for right to left. The default is to sweep from left to right first.
B. sweep_nbr: The total number of sweeps the user desires to make. If sweep_nbr= 0, then sweeps will continue until convergence.
C. norm: whether to normalize after each sweep or after all the sweeps. If norm = 0, the compressed MPS will be normalized after all the sweeps. If norm = 1, the compressed MPS will be normalized after each sweep.

* measurement.py defines the base measurement classe Measurement for all the measurement classes. It can compute the joint probability, correlation function, variance, and mean value at any time point. For the last two calculations, user can define specific values for being angry (up) or calm (down), and by default, they are 1 and -1 respectively. Both time and site number counts from 0, and a python-way of specifying the number from the end of list is also accepted. Each measurement task must be specified in a particular way, though measurement time can be omitted, in which case it is taken to be the most recent time computed. The class uses the solver to initialize, and the measure functions take the same input. They have different requested tasks and each measurement function should check the task properties. It will raise an exception if it does not fit. 
  For example, ("correlation", 1, [1, 2, 3, 4, 5]) measures <S1*S2*S3*S4*S5> at t = 1, task_type_string "Correlation", "Mean", "Variance", "Proba", more to be added some examples:("Correlation", 1, [1, 2, 3, 4, 5]) ("Mean", 1, [1]) ("Variance", 1, [1])("Proba", [(1, "up"), (3, "down")])It contains the following functions:
      1. addMeasureTask(self, task): The user can add different tasks by this function. The task is a tuple that shows up in the form (task_type_string, time_point, [points of interest]). Note in each task tuple, the integer following the task name will always be taken as the time step for all measurement. If time_point is not given, then it is assumed to be the last time step in calculation. It can also be negative which then follows Python's convention.
      2. measureProba(self, task): This function implements the measurement of the probability (possibly joint probability) for an event or several events to be realized
      3. measureCorrelation(self, task, up=None, down=None): This implements the measurement of the n-point correlation function
      4. measureMean(self, task, up=None, down=None): It calculates the mean value of the model.
      5. measureVariance(self, task, up=None, down=None): Variance is a special case of correlation, with only one variable. Just find the correct type and call measureCorrelation.
      6. getTimeTask(self, task, task_temp): For a given task, find the time for the computation and the true tasks. If the second element in the task is an integer, then it is taken as the time. If time is not given in the task at the second position, then it is the last one. task_temp is the list that stores true tasks.
* mpsmeasurement.py: It is the derived measurement class for MPS. It has the following assumptions: for a mps, physical index 0 represents down and physical index 1 represents up. solver.results has elements time ordered and start from time 0. the sites start from 0. The first mps is of shape (2,1,n), and the last mps of shape (2,n,1).
* exactmeasurement.py: It is the derived measurement classes for exact solution. 


This project is still in progress and is subject to the GPL 2 licence. Please feel free to contact the authors if you have some questions or suggestions.

Authors:
---
* Bin Xu (binarybin)
* Peiqi Wang (pqwang1026)
* Liangsheng Zhang (phzls)
* Peter Gjeltema (PJ)
* Jun Xiong (xiongPU)
