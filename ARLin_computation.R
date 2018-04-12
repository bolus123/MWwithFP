library(parallel)

#head.addr <- '/home/yuhuiyao/Documents/Github/MWwithFP/F&P.R'
head.addr <- 'C:/Dropbox/F&P/F&P.R'

#source(head.addr)

cores <- detectCores() - 1

rp <- 1000

m <- c(20)
n <- c(5)
cc <- 4.5625

shift.seq <- c(0, 0.01, 0.05)
scale.seq <- c(0.5, 0.75, 1)

x.sim <- 10
y.sim <- 10000
method <- 'Fligner-Policello'

cl <- makeCluster(cores)

clusterExport(cl, c('method', 'm', 'n', 'cc', 'x.sim', 'y.sim', 'head.addr'))
clusterCall(cl, function() source(head.addr))

res1 <- matrix(NA, nrow = length(shift.seq) * length(scale.seq), ncol = 4)
res2 <- matrix(NA, nrow = length(shift.seq) * length(scale.seq), ncol = 4)
res3 <- matrix(NA, nrow = length(shift.seq) * length(scale.seq), ncol = 4)
res4 <- matrix(NA, nrow = length(shift.seq) * length(scale.seq), ncol = 4)

i <- 0

for (shift in shift.seq) {

	for (scale in scale.seq) {

		clusterExport(cl, c('shift', 'scale'))

		i <- i + 1

		res1[i, 1] <- shift
		res1[i, 2] <- scale
		res1[i, 3] <- cc
		res1[i, 4] <- mean(unlist(parLapply(
			cl, 
			1:rp,
			function(x) {
				U.ARLin(m = m, n = n, cc = cc, method = method, FP.stat.type = 'exact', 
					x.sim = x.sim, y.sim = y.sim, X.dist = 'norm', Y.dist = 'norm', X.pars = c(0, 1), Y.pars = c(0, 1), 
					X.shift = 0, X.scale = 1, Y.shift = shift, Y.scale = scale)
			}
		)))


#############################################################################################################################
		res2[i, 1] <- shift
		res2[i, 2] <- scale
		res2[i, 3] <- cc
		res2[i, 4] <- mean(unlist(parLapply(
			cl, 
			1:rp,
			function(x) {
				U.ARLin(m = m, n = n, cc = cc, method = method, FP.stat.type = 'exact', 
					x.sim = x.sim, y.sim = y.sim, X.dist = 'norm', Y.dist = 't', X.pars = c(0, 1), Y.pars = c(5, 1), 
					X.shift = 0, X.scale = 1, Y.shift = shift, Y.scale = scale)
			}
		)))

#############################################################################################################################
		res3[i, 1] <- shift
		res3[i, 2] <- scale
		res3[i, 3] <- cc
		res3[i, 4] <- mean(unlist(parLapply(
			cl, 
			1:rp,
			function(x) {
				U.ARLin(m = m, n = n, cc = cc, method = method, FP.stat.type = 'exact', 
					x.sim = x.sim, y.sim = y.sim, X.dist = 'norm', Y.dist = 'laplace', X.pars = c(0, 1), Y.pars = c(0, 1), 
					X.shift = 0, X.scale = 1, Y.shift = shift, Y.scale = scale)
			}
		)))

#############################################################################################################################
		res4[i, 1] <- shift
		res4[i, 2] <- scale
		res4[i, 3] <- cc
		res4[i, 4] <- mean(unlist(parLapply(
			cl, 
			1:rp,
			function(x) {
				U.ARLin(m = m, n = n, cc = cc, method = method, FP.stat.type = 'exact', 
					x.sim = x.sim, y.sim = y.sim, X.dist = 'norm', Y.dist = 'gamma', X.pars = c(0, 1), Y.pars = c(5, 5), 
					X.shift = 0, X.scale = 1, Y.shift = shift, Y.scale = scale)
			}
		)))				

	}
}

stopCluster(cl)

#debug(U.ARLin)
#debug(dU.modified.f)
#debug(simulate.data)
#U.ARLin(m = 20, n = 5, cc = 4.5625, method = 'Fligner-Policello', FP.stat.type = 'exact', 
#			x.sim = 10, y.sim = 10000, X.dist = 'norm', Y.dist = 'norm', X.pars = c(0, 1), Y.pars = c(0, 1), 
#		X.shift = 0, X.scale = 1, Y.shift = 0.01, Y.scale = 1)

#a <- dU.modified.f(20, 5, type = 'exact', x.sim = 20, y.sim = 10000, 
#		X.dist = 'norm', Y.dist = 'norm', X.pars = c(0, 1), Y.pars = c(0, 1), 
#		X.shift = 0, X.scale = 1, Y.shift = 0, Y.scale = 1) 