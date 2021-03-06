library(parallel)

#head.addr <- '/home/yuhuiyao/Documents/Github/MWwithFP/F&P.R'
head.addr <- 'C:/Dropbox/F&P/F&P.R'


cores <- detectCores() - 2

rp <- 1000
#m.seq <- 200
m.seq <- c(25, 30, 35, 40, 45, 75, 125, 150, 200)
#m.seq <- 25
#m.seq <- c(25, 30)
n.seq <- c(5, 10)
ARL0 <- 370
#ARL0 <- 500
x.sim <- 10
y.sim <- 10000
tol <- .Machine$double.eps^0.25
method <- 'Mann-Whitney'

cl <- makeCluster(cores)

clusterExport(cl, c('method', 'ARL0', 'x.sim', 'y.sim', 'tol', 'head.addr'))
clusterCall(cl, function() source(head.addr))

res <- matrix(NA, nrow = length(m.seq) * length(n.seq), ncol = 4)

i <- 0

for (m in m.seq) {

	for (n in n.seq) {

		clusterExport(cl, c('m', 'n'))

		i <- i + 1

		res[i, 1] <- m
		res[i, 2] <- n
		res[i, 3] <- ARL0
		res[i, 4] <- mean(unlist(parLapply(
						cl, 
						1:rp, 
						function(x) {
							U.Charting.constants(m, n, ARL0, method = method, FP.stat.type = 'exact', 
								interval = c(1, 5), x.sim = x.sim, y.sim = y.sim, tol = tol)

						}
					)))


	}

}

stopCluster(cl)

proc.time()


#U.Charting.constants(20, 5, 370, method = 'Mann-Whitney', FP.stat.type = 'exact', 
#								interval = c(1, 5))