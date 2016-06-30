###### simulate data
set.seed(355)
source("functions.R")

#### test definitions of dsd
n <- 10
z <- sort(rnorm(n, sd = 20))
p <- rpois(n, lambda = 4)
p <- p / sum(p)
q <- rpois(n, lambda = 4)
q <- q / sum(q)

h <- maxGrad(z, p, q)
all.equal(h$dsd, sum(h$h * (q-p)))
myFlow <- minFlow(z, p, q)
flow <- myFlow$flow
all.equal(h$dsd, sum(flow * abs(outer(z, z, "-"))))


#### multivariate version and selection gradients
n <- 31
m <- 5  # number of traits
Z <- matrix(rnorm(m*n, sd = 1), ncol = m)
W <- exp(abs(Z[, 1]) - abs(Z[, 2])) + abs(Z[, 3] - 2*Z[, 4]) + (Z[, 5] - 1) ^ 2  + rexp(n, rate = 1)
tmp <- distSelGrad(Z, W)
# an alternative formula is lm(w*n - 1 ~ -1 + H)$coefficients, and since colSums(H) == 0 this is the same as lm(w*n ~ H)$coefficients[-1]
H <- tmp$matGrad
w <- W / sum(W)
all.equal(unname(lm(w*n ~ H - 1)$coefficients), tmp$sgrad)

(i <- sample(m, 1))  # random trait
zz <- Z[, i]
pp <- rep(1/n, n)
hh <- maxGrad(zz, pp, w)
c(hh$dsd, sum(hh$h * (w-pp)), tmp$DSD[i])


#### flow composition and decomposition
deComp <- flowDecomposition(myFlow)
D <- deComp$directional
N <- deComp$nondirectional
all.equal(D, flowDecomposition(list(flow = D, points = z))$directional)
all.equal(N, flowDecomposition(list(flow = N, points = z))$nondirectional)

all.equal(D + N - diag(colSums(D)), flow)
all.equal(D + N - diag(rowSums(N)), flow)
all.equal(rowSums(D), p)
all.equal(colSums(D), rowSums(N))
all.equal(colSums(N), q)

sum(D * abs(outer(z, z, "-")))
-sum(D * outer(z, z, "-"))		#### should be equal in absolute value
(s <- sum(z * (q-p)))			#### should be equal
all.equal(0, sum(N * outer(z, z, "-")))
sum(N * abs(outer(z, z, "-")))
sum(h$h * (q-p)) - abs(s)			#### should be equal



#### when pure directional selection takes place things there is a simple relationship between selection gradients and distributional selection gradients
D <- 5
s <- c(4, 3, -2, -5, 2)
P <- rWishart(1, D, diag(D))[, , 1] / D
beta <- solve(P) %*% s

pos <- 2 * (s >= 0) - 1   #### plus or minus one
d <- abs(s)
H <- P * outer(pos, pos, "*")
delta <- solve(H) %*% d
cbind(beta, delta)
