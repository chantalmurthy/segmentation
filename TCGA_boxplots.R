#http://stevencarlislewalker.wordpress.com/2012/06/28/overall-axis-labels-with-mfrow/
data <- read.table("C:/Users/Nameeta/Chantal/data/batch_protein_summaries_boxplots.txt", header=FALSE)

#average the repeated columns
HMOX1 <- as.numeric(as.character(data[-1,6])) + as.numeric(as.character(data[-1,7]))/2
TGFBI <- as.numeric(as.character(data[-1,4])) + as.numeric(as.character(data[-1,5]))/2
VCAM1 <- as.numeric(as.character(data[-1,2])) + as.numeric(as.character(data[-1,3]))/2
CD44 <- as.numeric(as.character(data[-1,8])) + as.numeric(as.character(data[-1,9]))/2

color_vec <- c("green", "green", "red", "red", "green")
labs <- c("Hlthy_src11", "Hlthy_src12", "GBM_src1", "GBM_src2", "Hlthy_src2")
par(mfrow=c(1,4), oma=c(4,4,2,0))
boxplot(HMOX1[which(data[,1]==labs[1])], HMOX1[which(data[,1]==labs[2])], HMOX1[which(data[,1]==labs[3])], HMOX1[which(data[,1]==labs[4])], HMOX1[which(data[,1]==labs[5])], range=0, col=color_vec, names=labs, notch=TRUE, outline=FALSE, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, boxwex=0.5, at=c(0,0.7,1.4,3,3.7),las=2)
stripchart(HMOX1[which(data[,1]==labs[1])], vertical=TRUE, add=TRUE, at=0, cex=1.6)
stripchart(HMOX1[which(data[,1]==labs[2])], vertical=TRUE, add=TRUE, at=0.7, cex=1.6)
stripchart(HMOX1[which(data[,1]==labs[3])], vertical=TRUE, add=TRUE, at=1.4, cex=1.6)
stripchart(HMOX1[which(data[,1]==labs[4])], vertical=TRUE, add=TRUE, at=3, cex=1.6)
stripchart(HMOX1[which(data[,1]==labs[5])], vertical=TRUE, add=TRUE, at=3.7, cex=1.6)
title("HMOX1", cex.main = 2)

boxplot(TGFBI[which(data[,1]==labs[1])], TGFBI[which(data[,1]==labs[2])], TGFBI[which(data[,1]==labs[3])], TGFBI[which(data[,1]==labs[4])], TGFBI[which(data[,1]==labs[5])], range=0, col=color_vec, names=labs, notch=TRUE, outline=FALSE, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, boxwex=0.5, at=c(0,0.7,1.4,3,3.7),las=2)
stripchart(TGFBI[which(data[,1]==labs[1])], vertical=TRUE, add=TRUE, at=0, cex=1.6)
stripchart(TGFBI[which(data[,1]==labs[2])], vertical=TRUE, add=TRUE, at=0.7, cex=1.6)
stripchart(TGFBI[which(data[,1]==labs[3])], vertical=TRUE, add=TRUE, at=1.4, cex=1.6)
stripchart(TGFBI[which(data[,1]==labs[4])], vertical=TRUE, add=TRUE, at=3, cex=1.6)
stripchart(TGFBI[which(data[,1]==labs[5])], vertical=TRUE, add=TRUE, at=3.7, cex=1.6)
title("TGFBI", cex.main = 2)

boxplot(VCAM1[which(data[,1]==labs[1])], VCAM1[which(data[,1]==labs[2])], VCAM1[which(data[,1]==labs[3])], VCAM1[which(data[,1]==labs[4])], VCAM1[which(data[,1]==labs[5])], range=0, col=color_vec, names=labs, notch=TRUE, outline=FALSE, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, boxwex=0.5, at=c(0,0.7,1.4,3,3.7),las=2)
stripchart(VCAM1[which(data[,1]==labs[1])], vertical=TRUE, add=TRUE, at=0, cex=1.6)
stripchart(VCAM1[which(data[,1]==labs[2])], vertical=TRUE, add=TRUE, at=0.7, cex=1.6)
stripchart(VCAM1[which(data[,1]==labs[3])], vertical=TRUE, add=TRUE, at=1.4, cex=1.6)
stripchart(VCAM1[which(data[,1]==labs[4])], vertical=TRUE, add=TRUE, at=3, cex=1.6)
stripchart(VCAM1[which(data[,1]==labs[5])], vertical=TRUE, add=TRUE, at=3.7, cex=1.6)
title("VCAM1", cex.main = 2)

boxplot(CD44[which(data[,1]==labs[1])], CD44[which(data[,1]==labs[2])], CD44[which(data[,1]==labs[3])], CD44[which(data[,1]==labs[4])], CD44[which(data[,1]==labs[5])], range=0, col=color_vec, names=labs, notch=TRUE, outline=FALSE, cex.axis=1.5, cex.names=1.5, cex.lab=1.5, boxwex=0.5, at=c(0,0.7,1.4,3,3.7),las=2)
stripchart(CD44[which(data[,1]==labs[1])], vertical=TRUE, add=TRUE, at=0, cex=1.6)
stripchart(CD44[which(data[,1]==labs[2])], vertical=TRUE, add=TRUE, at=0.7, cex=1.6)
stripchart(CD44[which(data[,1]==labs[3])], vertical=TRUE, add=TRUE, at=1.4, cex=1.6)
stripchart(CD44[which(data[,1]==labs[4])], vertical=TRUE, add=TRUE, at=3, cex=1.6)
stripchart(CD44[which(data[,1]==labs[5])], vertical=TRUE, add=TRUE, at=3.7, cex=1.6)
title("CD44", cex.main = 2)
title("Batch-based Protein Summaries", cex.main=3, outer=TRUE)

mtext('Z score ng/ml', side = 2, outer = TRUE, line = 2, cex.main=2)

