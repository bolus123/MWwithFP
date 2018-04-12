# Load the R MPI package

source(system.file("Rprofile", package = 'Rmpi'))

#loading functions

head.addr <- 'C:/Dropbox/F&P/F&P.R'

source(file = head.addr)

# loading settings

rp <- 10
m.seq <- c(20, 25)
#m.seq <- c(25, 30, 35, 40, 45, 75, 125, 150, 200)
n.seq <- c(5, 10)
ARL0.seq <- c(370, 500)
x.sim <- 10
y.sim <- 10000
tol <- .Machine$double.eps^0.25
method <- 'Fligner-Policello'

# Spawn N-1 workers

mpi.spawn.Rslaves(nslaves=mpi.universe.size()-1)

# Sending functions and settings to salves

mpi.bcast.Rfun2slave()

mpi.bcast.Robj2slave(x.sim)
mpi.bcast.Robj2slave(y.sim)
mpi.bcast.Robj2slave(tol)
mpi.bcast.Robj2slave(method)

# Set timer

start.time <- Sys.time()

res <- matrix(NA, nrow = length(m.seq) * length(n.seq) * length(ARL0.seq), ncol = 5)

i <- 0

# Computing CARL
for (m in m.seq) {

	for (n in n.seq) {

		for (ARL0 in ARL0.seq) {

			mpi.bcast.Robj2slave(m)
			mpi.bcast.Robj2slave(n)
			mpi.bcast.Robj2slave(ARL0)

			res[i, 1] <- m
			res[i, 2] <- n
			res[i, 3] <- ARL0
			res[i, 4] <- rp
			res[i, 5] <- mean(unlist(mpi.parLapply(
						1:rp
						,function(x) {
								U.Charting.constants(m, n, ARL0, method = method, FP.stat.type = 'exact', 
											interval = c(1, 15), x.sim = x.sim, y.sim = y.sim, tol = tol)
						}
					)))

		}
	}

}

end.time <- Sys.time()

res <- list(time = end.time - start.time, data = res.c.ii)

for (ii in 1:10000){

	add <- paste('/home/yyao17/2017Fall/CC/MPI/res.c.ii.MPI.', ii, '.Rdata', sep = '')

	if (!file.exists(add)) {
	    save(res, file = add)
	    break
	}


}

# Summarize CARL

#final.res.c.ii <- rep(0, length(unlist(res.c.ii[[1]])))

#for (i in 1:rp){

#	final.res.c.ii <- final.res.c.ii + res.c.ii[[i]]

#}

#final.res.c.ii <- final.res.c.ii / rp

#final.res.c.ii <- cbind(pars, final.res.c.ii)

#save(final.res.c.ii, file = '/home/yyao17/2017Fall/CC/MPI/final.res.c.ii.MPI.01112018.1.Rdata')

# Stop the worker processes

mpi.close.Rslaves()

# Close down the MPI processes and quit R

mpi.quit()


