##Makes plots of the projected values onto the first two principle components with each ROI colored differently
path <- "C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/"
feature_files <- list.files(path)
scalev <- 45
library('scatterplot3d')

for (file_ind in 1:length(feature_files)) {
	fileName <- feature_files[file_ind]
	sliceName0 <- strsplit(fileName, ".txt")
	sliceName <- strsplit(sliceName0[[1]], "feat_")[[1]][2]
	feature_matrix <- read.table(paste(path, fileName, sep=""))
	annotated_matrix <- feature_matrix[which(feature_matrix[,1] < 5),]
	annotations <- annotated_matrix[,1]
	#add normalization before doing PCA
	pca_input <- annotated_matrix[,-1]
	feature_pca <- prcomp(pca_input, scale = TRUE)
	colors <- matrix(data=0, nrow=nrow(pca_input), ncol=1)
	colors[which(annotations==1),1] <- "red" #CT
	colors[which(annotations==2),1] <- "green" #NE
	colors[which(annotations==3),1] <- "pink" #PAN
	colors[which(annotations==4),1] <- "purple"	#HBV
	png(paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/PCA/", sliceName, ".png", sep=""), height=700, width=1350)	
	plot(predict(feature_pca)[,1],predict(feature_pca)[,3],pch=21,col=colors, xlab="PC1", ylab="PC2")
	dev.off()
	
	# plot all companies loadings on the first, second, and third principal components and highlight points according to the sector they belong
	#s3d <- scatterplot3d(predict(feature_pca)[,1], predict(feature_pca)[,2], predict(feature_pca)[,3], xlab='Comp.1', ylab='Comp.2', zlab='Comp.3', color=colors, pch = 20)
	
	###MDS thing
	dat <- pca_input
	distmat = sqrt(1-cor(dat,use="na.or.complete")^2);   #  Or a different distance metric.  “dat” is the expression data matrix, typically in log2 space
	loc     = cmdscale(distmat,k=2, eig=TRUE);
	x       = -loc[[1]][,1];
	y       = loc[[1]][,2]
	loc1    = cmdscale(distmat,k=1, eig=TRUE);   # For percent variance explained calculation
	plot(x, y, pch=19, main=main,xlim=range(x)*1.2,ylim=range(y)*1.2,
	   ylab = paste("% Var explained, Dim. 2:  ",100*signif(loc$GOF[1]-loc1$GOF[1],3)),
	   xlab = paste("% Var explained, Dim. 1:  ",100*signif(loc1$GOF[1],3)))
	   
	#make plot of binary classification accuracy vs. number of PCs used to project the data
	
	#lda()
	
	
}
