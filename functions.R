maxGrad <- function(z, p, q)
####### calculate optimal h and the dsd
{
	if((n <- length(z)) != length(p)) stop("length of z and q need to be the same")
	if(n != length(q)) stop("length of z and p need to be the same")
	if(any(p < 0))	stop("p should be nonnegative")
	if(any(q < 0))	stop("q should be nonnegative")
	if(!isTRUE(all.equal(sp <- sum(p), sum(q))))	stop("p and q should have the same sum")
	if(!isTRUE(all.equal(sp, 1)))					warning("p does not sum up to one")
	
	if(UNSORTED <- is.unsorted(z))
	{
		or <- order(z)
		z <- z[or]
		p <- p[or]
		q <- q[or]
	}

	cumdiff <- (cumsum(p) - cumsum(q) )[-n]
	UNIQUE <- !any(0 == cumdiff)
	z.diff <- diff(z)
	hdiff <- z.diff * sign(cumdiff) # represent h(z[k]) - h(z[k-1])
	maxGrad <- cumsum(c(0, hdiff)) # represent h(z[k]);  arbitrary value for h(z[1])
	maxGrad <- maxGrad - sum(maxGrad * p) # make h have zero mean (with respect to the measure that gives pk mass to zk)
	
#	if(!isTRUE(all.equal(sum(z.diff * abs(cumdiff)  ), sum(maxGrad * (q - p)))))	stop("maximiser did not attain maximal value")

	list(h = if(UNSORTED)	maxGrad[order(or)]	else	maxGrad, dsd = sum(z.diff * abs(cumdiff)), UNIQUE = UNIQUE)
}

distSelGrad <- function(Z, W)
  ####### Z matrix of m traits of n individuals, W vector of fitness
{
  m <- dim(Z)
  if((n <- length(W)) != m[1]) stop("nrow(Z) should equal length(W)")
  if(any(W < 0))	stop("W should be nonnegative")
  if(all(W == 0))	stop("W should have at least one positive value")
  w <- W / sum(W)	#make it sum up to one (NOT have mean 1 like in text!)
  p <- rep(1/n, n)
  m <- m[2]	# number of traits
  DSD <- numeric(m)
  H <- Z		# maximisers
  for(i in 1:m)
  {
    res <- maxGrad(Z[, i], p = p, q = w)
    H[, i] <- res$h
    DSD[i] <- res$dsd
  }
  
  # the selection gradients are given by the below formula because colSums(H) == 0 so t(H) %*% H is proportional to cov(H)
  # an alternative formula is lm(w*n - 1 ~ -1 + H)$coefficients, and since colSums(H) == 0 this is the same as lm(w*n ~ H)$coefficients[-1]
  list(DSD = DSD, matGrad = H, sgrad = as.vector(n * solve(t(H) %*% H) %*% DSD))
}

minFlow <- function(z, p, q)
{
	if((n <- length(z)) != length(p)) stop("length of z and q need to be the same")
	if(n != length(q)) stop("length of z and p need to be the same")
	if(any(p < 0))	stop("p should be nonnegative")
	if(any(q < 0))	stop("q should be nonnegative")
	if(!isTRUE(all.equal(sp <- sum(p), sum(q))))	stop("p and q should have the same sum")
	if(!isTRUE(all.equal(sp, 1)))					warning("p does not sum up to one")

#### recast the optimal transportation problem as a linear program
#### this can be done using a more dedicated algorithm like the Hungarian method
	constraints <- matrix(0, nrow = 2*n, ncol = n^2)
	for(i in 1:n)
	{
		constraints[i, n*0:(n-1) + i] <- 1
		constraints[i + n, 1:n + (i-1)*n] <- 1
	}

#### remove last redundant constraint
	constraints <- constraints[-2*n, ]
	nconst <- 2*n - 1 # =nrow(constraints)


	cost.vec <- abs(outer(z, z, "-"))  #### L1 distance

	res <- lpSolve::lp(direction = "min", objective.in = cost.vec, const.mat = constraints, const.dir = rep("==", nconst), const.rhs = c(p, q[-n]))

	flow <- matrix(res$solution, ncol = n, nrow = n, byrow = FALSE)

#	if(!isTRUE(all.equal(res$objval, sum(diff(sort(z)) * abs(cumsum(p-q)[-n])  ))))	stop("minimiser did not attain minimal value")
	
	list(flow = flow, points = z)
}

flowDecomposition <- function(flow)
{
	z <- flow$points
	n <- length(z)
	flow <- flow$flow
	if(!isTRUE(all.equal(c(n, n), dim(flow))))	stop("flow$flow must be a square matrix of the same dimension as flow$points")
	if(any(0 > flow))	stop("flow must be nonnegative")

	if(is.unsorted(z))	stop("flow points must be sorted")#### adapt the code to this situation
	p <- rowSums(flow)
	q <- colSums(flow)
	L <- ! (H <- row(flow) < col(flow))
	l <- sum(   (flow * outer(z, z, "-"))[L]    )
	s <- sum(z * (q - p))
	if(s <= 0)
	{
		if(identical(l, 0))	D <- 0*flow else
		{
			
			D <- -flow*s/l
			D[H] <- 0
		}
	} else
	{
		h <- l + s
		D <- flow*s/h
		D[L] <- 0
	}

	D <- D + diag(p - rowSums(D))
	N <- flow - D
	N <- N + diag(q - colSums(N))
	
	list(directional = D, nondirectional = N)
}

composition <- function(D, N)
{
	if(!is.matrix(N) || !is.matrix(D) || !is.numeric(D) || !is.numeric(N))  stop("D and N should be numeric matrices")
	if(any(N < 0) || any(D < 0))	stop("D and N should be nonngative")
	n <- dim(N)
	if(n[1] != n[2])  stop("N needs to be a square matrix")
	if(!isTRUE(all.equal(n, dim(D))))  stop("D and N need to have the same dimension")
	n <- n[1]
	if(!isTRUE(all.equal(colSums(D), rowSums(N))))  stop("colSums(D) and rowSums(N) are different")
	if(   !isTRUE(all.equal(D[lower.tri(D)], rep(0, n*(n-1)/2)))
	   && !isTRUE(all.equal(D[upper.tri(D)], rep(0, n*(n-1)/2)))  )  warning("D is not triangular")
	
	D + N - diag(colSums(D))
}
