library(Hmisc)
#library(geoR)
#library(MASS)
#library(zoo)

library(parallel)

source('/home/yuhuiyao/Documents/Github/MWwithFP/F&P.R')

#boxcoxTrans <- function(x, lam1, lam2 = NULL) {
#
#    # if we set lambda2 to zero, it becomes the one parameter transformation
#    lam2 <- ifelse(is.null(lam2), 0, lam2)
#
#    if (lam1 == 0L) {
#      log(x + lam2)
#    } else {
#      (((x + lam2)^lam1) - 1) / lam1
#    }
# }

data1 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/LAB06HM_1999_2000.XPT")

data2 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/L06HM_B_2001_2002.XPT")

data3 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/L06HM_C_2003_2004.XPT")

data4 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/UHM_D_2005_2006.XPT")

data5 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/UHM_E_2007_2008.XPT")

data6 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/UHM_F_2009_2010.XPT")

data7 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/UHM_G_2011_2012.XPT")

data8 <- sasxport.get("/home/yuhuiyao/Documents/F&P/UrineData/UM_H_2013_2014.XPT")



ref.samp <- data1$urxupb[is.na(data1$urxupb) == FALSE]

test.samp1 <- data2$urxupb[is.na(data2$urxupb) == FALSE]
test.samp2 <- data3$urxupb[is.na(data3$urxupb) == FALSE]
test.samp3 <- data4$urxupb[is.na(data4$urxupb) == FALSE]
test.samp4 <- data5$urxupb[is.na(data5$urxupb) == FALSE]
test.samp5 <- data6$urxupb[is.na(data6$urxupb) == FALSE]
test.samp6 <- data7$urxupb[is.na(data7$urxupb) == FALSE]
test.samp7 <- data8$urxupb[is.na(data8$urxupb) == FALSE]

U1 <- U.stat.f(ref.samp, test.samp1)
U2 <- U.stat.f(ref.samp, test.samp2)
U3 <- U.stat.f(ref.samp, test.samp3)
U4 <- U.stat.f(ref.samp, test.samp4)
U5 <- U.stat.f(ref.samp, test.samp5)
U6 <- U.stat.f(ref.samp, test.samp6)
U7 <- U.stat.f(ref.samp, test.samp7)

U.modified1 <- U.modified.stat.f(ref.samp, test.samp1, type = 'exact') 
U.modified2 <- U.modified.stat.f(ref.samp, test.samp2, type = 'exact') 
U.modified3 <- U.modified.stat.f(ref.samp, test.samp3, type = 'exact') 
U.modified4 <- U.modified.stat.f(ref.samp, test.samp4, type = 'exact') 
U.modified5 <- U.modified.stat.f(ref.samp, test.samp5, type = 'exact') 
U.modified6 <- U.modified.stat.f(ref.samp, test.samp6, type = 'exact') 
U.modified7 <- U.modified.stat.f(ref.samp, test.samp7, type = 'exact') 


cc.FP1 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp1), ARL0 = 370)
cc.FP2 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp2), ARL0 = 370)
cc.FP3 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp3), ARL0 = 370)
cc.FP4 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp4), ARL0 = 370)
cc.FP5 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp5), ARL0 = 370)
cc.FP6 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp6), ARL0 = 370)
cc.FP7 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp7), ARL0 = 370)

cc.MW1 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp1), ARL0 = 370, method = 'Mann-Whitney')
cc.MW2 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp2), ARL0 = 370, method = 'Mann-Whitney')
cc.MW3 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp3), ARL0 = 370, method = 'Mann-Whitney')
cc.MW4 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp4), ARL0 = 370, method = 'Mann-Whitney')
cc.MW5 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp5), ARL0 = 370, method = 'Mann-Whitney')
cc.MW6 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp6), ARL0 = 370, method = 'Mann-Whitney')
cc.MW7 <- U.Charting.constants(m = length(ref.samp), n = length(test.samp7), ARL0 = 370, method = 'Mann-Whitney')



#pars.ref.samp <- boxcoxfit(ref.samp)
#
#trans.ref.samp <- boxcoxTrans(ref.samp, lam1 = as.numeric(pars.ref.samp$lambda))
#
#p1 <- hist(trans.ref.samp, breaks=20)
#
#breaks_cdf <- pnorm(p1$breaks, pars.ref.samp$beta.normal, 
#	sqrt(pars.ref.samp$sigmasq.normal))
#null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
#
#a <- chisq.test(p1$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)