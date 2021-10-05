library(difR)
library(ltm)
library(mirt)
setwd("C:/Users/sjoo/Dropbox/Research/PPMC DIF/Real Data")

data1 <- read.table("Data1.dat", quote="\"", comment.char="")
data2 <- read.table("Data2.dat", quote="\"", comment.char="")

ref_x <- cbind(rep(1,nrow(data1)),data1)
foc_x <- cbind(rep(2,nrow(data2)),data2)
colnames(ref_x) <- colnames(foc_x) <- c("grp",paste("I",seq(1:10),sep=""))

fit1 <- mirt(ref_x[,2:11],1,itemtype="2PL",)
fit2 <- mirt(foc_x[,2:11],1,itemtype="2PL")
M2(fit1)
M2(fit2)

x <- rbind(ref_x,foc_x)
grp <- as.character(x[,1])
multgroup <- multipleGroup(x[,2:11], 1, group = grp, invariance=c('slopes', 'intercepts'))
coef(multgroup, simplify=TRUE, IRTpars = TRUE)
      
# Lord's chisquare test
lordtest <- difLord(x[,2:11], group = x[,"grp"], focal.name = 2, model = "2PL")
# Raju's DIF test
rajutest <- difRaju(x[,2:11], group = x[,"grp"], focal.name = 2, model = "2PL")

lordtest
rajutest


