# Distributional selection differentials

R code for computations related to distributional selection differentials.

Requires the R package lpSolve for the calculation of the optimal flow.

## Contents
### The file functions.R contains the following functions.

| Function | Description |
| --- | --- |
| maxGrad | calculates the (standardised) maximiser h  |
| minFlow | calculates the optimal flow |
| flowDecomposition | calculates the decomposition of a given flow into a directional and nondirectional flow |
| composition | inverse function of flowDecomposition |

### tests.R
This file contains some simulations for trying out the functions in the file functions.R

### Copyright (c) 2016 Yoav Zemel and Jonathan M. Henshaw
