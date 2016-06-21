###### simulate data
#set.seed(355)
n <- 10
z <- sort(rnorm(n, sd = 20))
p <- rpois(n, lambda = 4)
p <- p / sum(p)
q <- rpois(n, lambda = 4)
q <- q / sum(q)

h <- maxGrad(z, p, q)
myFlow <- minFlow(z, p, q)
flow <- myFlow$flow
deComp <- flowDecomposition(myFlow)
D <- deComp$directional
N <- deComp$nondirectional


#### some tests
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
sum(h * (q-p)) - abs(s)			#### should be equal



D <- 5
s <- c(4, 3, -2, -5, 2)
P <- rWishart(1, D, diag(D))[, , 1] / D
beta <- solve(P) %*% s

pos <- 2 * (s >= 0) - 1   #### plus or minus one
d <- abs(s)
H <- P * outer(pos, pos, "*")
delta <- solve(H) %*% d
cbind(beta, delta)
