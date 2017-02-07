# Soot-Consumption-Modeling-with-Bayesian-Statistics
All the files are supplementary material to the paper that bears this repository's name.

The first file, the Data file, is found in two forms. The first form is a .xlsx file which contains all the data used in the paper's parameter calibration in a raw form along with derivations done to convert that raw data to the instumental model used in the soot consumption model. The three .csv files contain the same information but in a different format and in three files rather than one.

The second file, the bayesian.m file, is a Matlab function file which perform a Bayesian statistical evaluation. The input to this file are: a sample space of parameters to be evaluated, a prior multidimensional PDF, data against which the parameters will be evaluated, and the model to which these parameters belong. Note that when this function is used, it can be subject to truncation errors, due to very small likelihoods. It is recommended that one should not process more than 25 data points at once.
The Tutorial.m file contains a brief tutorial of how the bayesian.m function file is used. The tutorial evaluates some synthetic data and calibrate parameters in a quadratic model (found in the quadratic_model.m file) to fit that data.
