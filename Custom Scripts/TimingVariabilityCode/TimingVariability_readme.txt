TIMING VARIABILITY MODEL CODE

Timing variability model code associated with Glaze & Troyer (2012). This document gives an overview of the functions in this folder. See manuscript for model description and derivation, and comments in each function for details about argument structure. 


This folder contains 3 main functions and 2 auxillary functions:

MAIN FUNCTIONS
timing_var_EM.m: Runs EM algorithm on a single set of initial parameter estimates, with convergence criteria set either by user or by default. Returns parameter and latent variable estimates, as well as a few other pieces of model information.

timing_var.m: Runs timing_var_EM.m M times for M different initial parameter estimates, convergence criteria set either by user or default. Procedure for determining optimal parameter set is divided into 2 stages as described in manuscript.

timing_var_BIC.m: Runs timing_var.m while varying the dimensionality of the global variable, with the set of dimensions determined by the user. Determines optimal dimensionality using the Bayesian Information Criterion and returns a structure containing model information associated with that dimensionality (parameter estimates, etc.).


AUXILLARY FUNCTIONS
rotate_max.m: Rotates global weight matrix to maximize the sum of entries in the first column, orthogonalizes remaining columns (as described in manuscript).

args2timing_vars.m: Assigns default values to a give variable when user has not.