# dsd
Distributional selection differentials

R code for computations related to distributional selection differentials

functions.R:  maxGrad calculates the (standardised) maximiser h
              minFlow calculates the optimal flow, using an LP solver
              flowDecomposition calculates the decomposition of a given flow into a directional and nondirectional flow
              composition is the inverse of flowDecomposition

tests.R       some simulations for trying out the functions in the file functions.R
