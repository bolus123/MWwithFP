##############################################################################
#
##############################################################################

#library(NSM3)

##############################################################################

#U.stat.f <- function(X, Y){



#}

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi) 

W.stat.f <- function(X, Y) {

	m <- length(X)
	n <- length(Y)

	XY <- c(X, Y)

	id <- 1:(m + n)

	r <- rank(XY)

	Z <- id <= m

	sum(r[Z])

}

U.stat.f <- function(X, Y) {

	m <- length(X)
	n <- length(Y)

	U <- W.stat.f(X, Y) - m / 2 * (m + 1)

	U <- m * n - U

	(U  - m * n / 2) / sqrt((m * n * (m + n + 1) / 12))

}


dU.f <- function(m, n, x.sim = 20, y.sim = 10000) {

	X <- matrix(rnorm(m * x.sim), nrow = x.sim, ncol = m)

	dU <- matrix(NA, nrow = x.sim, ncol = y.sim)

	for (i in 1:x.sim) {
		for (j in 1:y.sim) {
			Y <- rnorm(n)
			dU[i, j] <- U.stat.f(X[i, ], Y)
		}
	}


	return(dU)
}

U.modified.stat.f <- function(X, Y, type = 'exact') {

	m <- length(X)
	n <- length(Y)

	rank.vec <- rank(c(X, Y))

    x.tmp <- rank.vec[1:m]
    y.tmp <- rank.vec[(m + 1):(m + n)]
    p.vec <- unlist(lapply(x.tmp, function(x) {
        sum(x > y.tmp) + 0.5 * sum(x == y.tmp)
    }))
    q.vec <- unlist(lapply(y.tmp, function(y) {
        sum(y > x.tmp) + 0.5 * sum(y == x.tmp)
    }))
    p.bar <- mean(p.vec)
    q.bar <- mean(q.vec)
    v1 <- sum((p.vec - p.bar)^2)
    v2 <- sum((q.vec - q.bar)^2)

    if (type == 'exact') {

    	res <- (sum(q.vec) - sum(p.vec))/(2 * sqrt((n - 1) / n * v1 + (m - 1) / m * v2 + 
        p.bar * q.bar))

    } else if (type == 'approx') {

    	res <- (sum(q.vec) - sum(p.vec))/(2 * sqrt(v1 + v2 + p.bar * q.bar))
    }

    return(res)

}

dU.modified.f <- function(m, n, type = 'exact', x.sim = 20, y.sim = 10000) {

	X <- matrix(rnorm(m * x.sim), nrow = x.sim, ncol = m)
	#Y <- matrix(rnorm(), )

	dU.modified <- matrix(NA, nrow = x.sim, ncol = y.sim)

	for (i in 1:x.sim) {
		for (j in 1:y.sim) {
			Y <- rnorm(n)
			dU.modified[i, j] <- U.modified.stat.f(X[i, ], Y, type = type) 
		}
	}

	return(dU.modified)

}

U.Charting.constants <- function(m, n, ARL0, method = 'Fligner-Policello', FP.stat.type = 'exact', interval = c(1, 3), x.sim = 20, y.sim = 10000, tol = .Machine$double.eps^0.25){

	root.finding <- function(cc, y.sim, ARL0, stat){

		ARLin <- mean( apply(stat, 1, function(x) sum(x < -cc| x > cc) / y.sim) ^ -1 )

		#cat('cc:', cc, 'ARLin', ARLin, '\n')

		ARL0 - ARLin

	}

	if (method == 'Mann-Whitney') {

		stat <- dU.f(m, n, x.sim = x.sim, y.sim = y.sim) 

	} else if (method == 'Fligner-Policello') {

		stat <- dU.modified.f(m, n, type = FP.stat.type, x.sim = x.sim, y.sim = y.sim) 

	}

	uniroot(root.finding, interval = interval, y.sim = y.sim, ARL0 = ARL0, stat = stat, tol = tol)$root

}

U.control.chart <- function(X, Y, cc = NA, method = 'Fligner-Policello', plot.option = TRUE, ARL0 = 370, 
	FP.stat.type = 'exact', interval = c(1, 3), x.sim = 20, y.sim = 10000, tol = .Machine$double.eps^0.25){

	m <- length(X)
	n <- dim(Y)[2]
	r <- dim(Y)[1]

	if (method == 'Fligner-Policello') {

		stat <- apply(Y, 1, function(z) U.modified.stat.f(X, z, type = FP.stat.type))

	} else if (method == 'Mann-Whitney') {

		stat <- apply(Y, 1, function(z) U.stat.f(X, z))

	}

	if (is.na(cc)) {

		UCL <- U.Charting.constants(m, n, ARL0, method = method, FP.stat.type = FP.stat.type, 
			interval = interval, x.sim = x.sim, y.sim = y.sim, tol = tol)
		LCL <- -UCL


	} else {

		UCL <- cc
		LCL <- -cc

	}

	
	if (plot.option == TRUE) {
		if (method == 'Fligner-Policello') {
			plot(c(1, r), c(min(stat, LCL), max(stat, UCL)), type = 'n', 
				xlab = 'Subgroups', ylab = 'Charting Statistics', 
				main = 'Mann-Whitney Control Chart with Fligner-Policello Replacement')
		} else if (method == 'Mann-Whitney') {
			plot(c(1, r), c(min(stat, LCL), max(stat, UCL)), type = 'n', 
				xlab = 'Subgroups', ylab = 'Charting Statistics', main = 'Mann-Whitney Control Chart')
		}
		points(1:r, stat, type = 'o', pch = 1)
		abline(h = UCL, lty = 2, col = 'red')
		abline(h = LCL, lty = 2, col = 'red')
		abline(h = 0, lty = 2)
	}

	res <- list(Charting.stat = stat, CL = 0, LCL = LCL, UCL = UCL)

	return(res)

}

Welch.t.f0 <- function(X, Y, ARL0 = 370, c4.option = TRUE, plot.option = TRUE) {

	m <- length(X)
	n <- dim(Y)[2]
	r <- dim(Y)[1]

	X.bar <- mean(X)
	Y.bar.vec <- rowMeans(Y)

	S2.X <- var(as.vector(X))
	S2.Y.vec <- diag(var(t(Y)))

	eta.hat.vec <- S2.Y.vec / S2.X

	S2.vec <- S2.X / m + eta.hat.vec * S2.X / n

	nu.numer <- (1 + eta.hat.vec * m / n) ^ 2
	nu.denom <- 1 / (m - 1) + eta.hat.vec ^ 2 * m ^ 2 / n ^ 2 / (n - 1)

	nu <- nu.numer / nu.denom

	if (c4.option == TRUE) {
		for (i in 1:r) {

			c4.vec <- rep(NA, r)	
			c4.vec[i] <- ifelse(c4.option == TRUE, c4.f(nu[i]), 1)

		}
		
	} else {
		c4.vec <- rep(1, r)
	}

	Charting.stat <- (Y.bar.vec - X.bar) / sqrt(S2.vec) * c4.vec

	quantile.t.vec <- qt(1 - 1 / 2 / ARL0, nu)

	UCL <- quantile.t.vec * c4.vec
	LCL <- -UCL

	res <- list(Charting.stat = Charting.stat, LCL = LCL, UCL = UCL)

	return(res)

}
#debug(Welch.t.f0)
#Welch.t.f0(X, Y)

Welch.t.f1 <- function(X, Y, ARL0 = 370, c4.option = TRUE, plot.option = TRUE) {

	m <- length(X)
	n <- dim(Y)[2]
	r <- dim(Y)[1]

	X.bar <- mean(X)
	Y.bar.vec <- rowMeans(Y)

	S2.X <- var(as.vector(X))
	S2.Y.vec <- diag(var(t(Y)))

	eta.hat.vec <- S2.Y.vec / S2.X

	nu.numer <- (1 + eta.hat.vec * m / n) ^ 2
	nu.denom <- 1 / (m - 1) + eta.hat.vec ^ 2 * m ^ 2 / n ^ 2 / (n - 1)

	nu <- nu.numer / nu.denom

	quantile.t.vec <- qt(1 - 1 / 2 / ARL0, nu)

	c4 <- ifelse(c4.option == TRUE, c4.f(m - 1), 1)

	factor.Kwt <- c4 * sqrt(1 + eta.hat.vec * m / n)

	Kwt <- quantile.t.vec * factor.Kwt

	Charting.stat <- (Y.bar.vec - X.bar) / sqrt(S2.X / m) * c4

	LCL <- -Kwt
	UCL <- Kwt

	if (plot.option == TRUE) {

		plot(c(1, r), c(min(Charting.stat, LCL), max(Charting.stat, UCL)), type = 'n', 
				xlab = 'Subgroups', ylab = 'Charting Statistics', 
				main = 'Welch t Control Chart')

		points(1:r, Charting.stat, type = 'o', pch = 1)
		points(1:r, LCL, lty = 2, col = 'red', type = 'l')
		points(1:r, UCL, lty = 2, col = 'red', type = 'l')
		abline(h = 0, lty = 2)


	}

	res <- list(Charting.stat = Charting.stat, CL = 0, LCL = LCL, UCL = UCL)

	return(res)

}

#debug(Welch.t.f)
#Welch.t.f1(X, Y)

Welch.t.charting.constants <- function(m, n, eta, ARL0 = 370, interval = c(1, 3), c4 = 1) {

	integrand.f <- function(g, m, n, eta, Kwt, c4) {

		nu.numer <- (1 + g * eta * m / n) ^ 2
		nu.denom <- 1 / (m - 1) + (g * eta) ^ 2 * m ^ 2 / n ^ 2 / (n - 1)

		nu <- nu.numer / nu.denom

		prob <- (pt(Kwt / c4 / sqrt(1 + g * eta * m / n), nu) 
			- pt(-Kwt / c4 / sqrt(1 + g * eta * m / n), nu)) * df(g, n - 1, m - 1)

		#ARL <- (1 - prob) ^ - 1 * df(g, n - 1, m - 1)

		#return(ARL)

		return(prob)

	}

	integral.f <- function(m, n, eta, Kwt, c4){

		integrate(integrand.f, lower = 0, upper = Inf, m = m, n = n, eta = eta, Kwt = Kwt, c4 = c4)$value

	}

	root.finding <- function(Kwt, ARL0, m, n, eta, c4){

		ARLin <- (1 - integral.f(m, n, eta, Kwt, c4)) ^ -1

		cat('ARLin', ARLin, '\n')

		ARL0 - ARLin

	}

	m <- length(X)
	n <- dim(Y)[2]
	r <- dim(Y)[1]	

	uniroot(root.finding, interval = interval, m = m, n = n, eta = eta, c4 = c4, ARL0 = ARL0)$root


}

#Welch.t.charting.constants(125, 5, 1, interval = c(5, 6))

Welch.t.f2 <- function(X, Y, eta = 1, ARL0 = 370, interval = c(1, 3), c4.option = TRUE, plot.option = TRUE) {

	m <- length(X)
	n <- dim(Y)[2]
	r <- dim(Y)[1]

	X.bar <- mean(X)
	Y.bar.vec <- rowMeans(Y)

	S2.X <- var(as.vector(X))

	c4 <- ifelse(c4.option == TRUE, c4.f(m - 1), 1)

	Kwt <- Welch.t.charting.constants(m, n, eta, ARL0 = ARL0, interval = interval, c4 = c4)

	Charting.stat <- (Y.bar.vec - X.bar) / sqrt(S2.X / m) * c4

	LCL <- -Kwt
	UCL <- Kwt

	if (plot.option == TRUE) {

		plot(c(1, r), c(min(Charting.stat, LCL), max(Charting.stat, UCL)), type = 'n', 
				xlab = 'Subgroups', ylab = 'Charting Statistics', 
				main = 'Welch t Control Chart')

		points(1:r, Charting.stat, type = 'o', pch = 1)
		abline(h = LCL, lty = 2, col = 'red')
		abline(h = UCL, lty = 2, col = 'red')
		abline(h = 0, lty = 2)


	}

	res <- list(Charting.stat = Charting.stat, CL = 0, LCL = LCL, UCL = UCL)

	return(res)

}

#Welch.t.f2(X, Y, eta = 1, ARL0 = 370, interval = c(1, 50), c4.option = TRUE, plot.option = TRUE) 


#rp <- 120
#
#cc.vec <- rep(NA, rp)
#
#for (i in 1:rp){
#
#	cc.vec[i] <- U.modified.charting.constants(125, 5, 370, interval = c(1, 10), x.sim = 100, y.sim = 10000)
#
#}


#X <- matrix(c(
#	1.3235,	1.4128,	1.6744,	1.4573,	1.6914,
#	1.4314,	1.3592,	1.6075,	1.4666,	1.6109,
#	1.4284,	1.4871,	1.4932,	1.4324,	1.5674,
#	1.5028,	1.6352,	1.3841,	1.2831,	1.5507,
#	1.5604,	1.2735,	1.5265,	1.4363,	1.6441,
#	1.5955,	1.5451,	1.3574,	1.3281,	1.4198,
#	1.6274,	1.5064,	1.8366,	1.4177,	1.5144,
#	1.419,	1.4303,	1.6637,	1.6067,	1.5519,
#	1.3884,	1.7277,	1.5355,	1.5176,	1.3688,
#	1.4039,	1.6697,	1.5089,	1.4627,	1.522,
#	1.4158,	1.7667,	1.4278,	1.5928,	1.4181,
#	1.5821,	1.3355,	1.5777,	1.3908,	1.7559,
#	1.2856,	1.4106,	1.4447,	1.6398,	1.1928,
#	1.4951,	1.4036,	1.5893,	1.6458,	1.4969,
#	1.3589,	1.2863,	1.5996,	1.2497,	1.5471,
#	1.5747,	1.5301,	1.5171,	1.1839,	1.8662,
#	1.368,	1.7269,	1.3957,	1.5014,	1.4449,
#	1.4163,	1.3864,	1.3057,	1.621,	1.5573,
#	1.5796,	1.4185,	1.6541,	1.5116,	1.7247,
#	1.7106,	1.4412,	1.2361,	1.382,	1.7601,
#	1.4371,	1.5051,	1.3485,	1.567,	1.488,
#	1.4738,	1.5936,	1.6583,	1.4973,	1.472,
#	1.5917,	1.4333,	1.5551,	1.5295,	1.6866,
#	1.6399,	1.5243,	1.5705,	1.5563,	1.553,
#	1.5797,	1.3663,	1.624,	1.3732,	1.6887
#), nrow = 25, ncol = 5, byrow = TRUE)
#
#Y <- matrix(c(
#	1.4483,	1.5458,	1.4538,	1.4303,	1.6206,
#	1.5435,	1.6899,	1.583,	1.3358,	1.4187,
#	1.5175,	1.3446,	1.4723,	1.6657,	1.6661,
#	1.5454,	1.0931,	1.4072,	1.5039,	1.5264,
#	1.4418,	1.5059,	1.5124,	1.462,	1.6263,
#	1.4301,	1.2725,	1.5945,	1.5397,	1.5252,
#	1.4981,	1.4506,	1.6174,	1.5837,	1.4962,
#	1.3009,	1.506,	1.6231,	1.5831,	1.6454,
#	1.4132,	1.4603,	1.5808,	1.7111,	1.7313,
#	1.3817,	1.3135,	1.4953,	1.4894,	1.4596,
#	1.5765,	1.7014,	1.4026,	1.2773,	1.4541,
#	1.4936,	1.4373,	1.5139,	1.4808,	1.5293,
#	1.5729,	1.6738,	1.5048,	1.5651,	1.7473,
#	1.8089,	1.5513,	1.825,	1.4389,	1.6558,
#	1.6236,	1.5393,	1.6738,	1.8698,	1.5036,
#	1.412,	1.7931,	1.7345,	1.6391,	1.7791,
#	1.7372,	1.5663,	1.491,	1.7809,	1.5504,
#	1.5971,	1.7394,	1.6832,	1.6677,	1.7974,
#	1.4295,	1.6536,	1.9134,	1.7272,	1.437,
#	1.6217,	1.822,	1.7915,	1.6744,	1.9404
#
#), nrow = 20, ncol = 5, byrow = TRUE)
#
#
#
##U.control.chart(X, Y, 7.5899)
#
#U.Charting.constants(25, 5, 500, method = 'Fligner-Policello', FP.stat.type = 'exact', 
#	interval = c(1, 10), x.sim = 10, y.sim = 10000, tol = .Machine$double.eps^0.25)
#
#
#U.control.chart(X, Y, cc = NA, method = 'Fligner-Policello', plot.option = TRUE, ARL0 = 500, 
#	FP.stat.type = 'exact', interval = c(1, 9), x.sim = 20, y.sim = 10000, tol = .Machine$double.eps^0.25)
#
#
#U.control.chart(X, Y, cc = NA, method = 'Mann-Whitney', plot.option = TRUE, ARL0 = 370, 
#	interval = c(1, 3), x.sim = 20, y.sim = 10000, tol = .Machine$double.eps^0.25)#