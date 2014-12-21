MPS on Probability
===

This program studies the stochastic process problems with an algorithm called the matrix product state (MPS) that showed its power first in many body quantum physics. Noticing the analogy between probability theory and quantum mechanics, we adapted the algorithm to deal with many variable stochastic process problems.

We can use transitional matrix methods to study stochastic processes and extract a lot of interesting properties. However, just like in the many body quantum system, the transitional matrix grows exponentially in a system of multiple degrees of freedom and becomes intractible with 30~50 variables. However, many models benefit from the property of locality and it is proved that, writing states as product of small matrices can effectively reduce the system complexity. We write transitional matrices as well as state vectors in terms of product of matrices and solve the time evolution problem efficiently.

A detailed description of the problems, models and the algorithm is in the documentation file Project_note.tex and a detailed user's manual will be available before the mid-January, 2015. If you want to play with the program before the beta version is released. 

Here is a quick start guide.
---
proba.py is the main caller and you can tune the parameters in the main function
model.py defines some models 
exactsolver.py and mpssolver.py define two solvers: the former is the traditional transitional matrix method and the latter is the MPS solver. You can compare them for small sizes.
exactmeasurement.py and mpsmeasurement.py are measurement classes that can compute mean, variance, correlation functions and joint probability.

This project is still in progress and is subject to the GPL 2 licence. Please feel free to contact the authors if you have some questions or suggestions.

Authors:
Bin Xu (binarybin)
Peiqi Wang (pqwang1026)
Liangsheng Zhang (phzls)
Peter Gjeltema (PJ)
Jun Xiong (xiongPU)
