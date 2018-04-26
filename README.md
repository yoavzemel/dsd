# Distributional selection differentials

## IMPORTANT:  THIS VERSION OF THE CODE IS OLD.  PLEASE USE JJONO'S UPDATED VERSION AT https://github.com/jjono/DSD

R code for computations related to distributional selection differentials, defined in the paper:

Henshaw JM, Zemel Y (2017). A unified measure of linear and nonlinear selection on quantitative traits. Methods in Ecology and Evolution 8(5): 604-614 (doi:10.1111/2041‑210X.12685)

Requires the R package lpSolve for the calculation of the optimal flow.


## Contents
### The file functions.R contains the following functions.

| Function | Description |
| --- | --- |
| distSelGrad | calculates the distributional selection gradients and the distributional selection differentials for a given matrix of trait values and a given vector of (absolute) fitness  |
| maxGrad | calculates the (standardised) maximiser h  |
| minFlow | calculates the optimal flow |
| flowDecomposition | calculates the decomposition of a given flow into a directional and nondirectional flow |
| composition | inverse function of flowDecomposition |

### tests.R
This file contains some simulations for trying out the functions in the file functions.R

### Copyright (c) 2016 Yoav Zemel and Jonathan M. Henshaw
