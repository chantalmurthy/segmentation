raw_dat <- read.table("C:/Users/Nameeta/Chantal/results/TCGA_June2014/normalized_KIRP.txt", header=TRUE)
genes <- raw_dat[,1]
dat <- raw_dat[,-1]
mads=apply(dat,1,mad)
dat=dat[rev(order(mads))[1:5000],]
num_samples <- ncol(dat)
dat = sweep(dat,1, apply(dat,1,median,na.rm=T))
library(ConsensusClusterPlus)
title=tempdir()
maxK = 6
reps=50
pItem=0.8
pFeature=1
results = ConsensusClusterPlus(as.matrix(dat),maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature, title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png") #what is the difference between matrix, data frame and data set?

#later understand how "cluster consensus/item consensus" computed
#icl = calcICL(results,title=title,plot="png")
#icl[["clusterConsensus"]]

#visualize result

#Plot Consensus Matrix CDFs for different K
plot(0,0,xlim=c(0,1),ylim=c(0,1),type="o")
range <- seq(0.1,1,0.1)
colors <- rainbow(maxK-1)
CDF_mat <- matrix(data=0, nrow=maxK-1, ncol=length(range))
dif_area_under_CDF <- c()
for (k in 2:maxK) {
	mat <- results[[k]]['consensusMatrix']$consensusMatrix
	area_under_CDF <- 0
	for (ind in 1:length(range)) {
		val <- length(which(mat < range[ind]))/(nrow(mat)*ncol(mat))
		CDF_mat[k-1,ind] <- val
		if (ind > 1) {
			area_under_CDF <- area_under_CDF + 0.1*(val - CDF_mat[k-1,ind-1])
		}
	}
	points(range, CDF_mat[k-1,], col=colors[k-1],type="o")
	dif_area_under_CDF <- c(dif_area_under_CDF, area_under_CDF)
}

#difference between change in area under CDF for each pair of consecutive cluster numbers between 2 to maxK-1
plot(c(2:maxK), dif_area_under_CDF,type="o")

#other plots



