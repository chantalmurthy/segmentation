directory <- "C:/Users/Nameeta/Documents/PR002/sorted_files/"
file_names <- list.files(directory)
map <- new.env(hash=T, parent=emptyenv())
map2 <- new.env(hash=T, parent=emptyenv())
patient_names <- c()
#patient name -block -section
for (i in 1:length(file_names)) {
	filename <- file_names[i]
	splt <- strsplit(filename, "_")
	patient <- splt[[1]][1]
	patient2 <- strsplit(patient, "-")[[1]]
	patient_name <- paste(patient2[1], patient2[2], patient2[3], sep="-")
	patient_names <- c(patient_names, patient_name)
	temp <- strsplit(patient, "\\.")[[1]][1]
	block_name <- substr(temp, nchar(temp), nchar(temp))
	if (is.null(map[[patient_name]])) {
		map[[patient_name]] <- new.env(hash=T, parent=emptyenv())
	}
	map[[patient_name]][[block_name]] <- c(map[[patient_name]][[block_name]], filename)
	map2[[patient_name]] <- c(map2[[patient_name]], filename)
}

block_names <- c(LETTERS, "AA", "AB", "AC", "AD")
for (i in 1 : length(patient_names)) {
	print(paste("patient: ", patient_names[i], sep=""))
	#make directory for patient
	patient <- patient_names[i]
	dir.create(file.path(paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, sep="")), showWarnings = FALSE)
	dir.create(file.path(paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations", sep="")), showWarnings = FALSE)
	dir.create(file.path(paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/RelNucFraction", sep="")), showWarnings = FALSE)
	num_annotations_patient <- matrix(data=0, nrow=1, ncol=7) #CT, NE, HBV, MVP, LE, IT, PAN
	num_annotations_z <- matrix(data=0, nrow=40, ncol=7) #Just guessing that number of z levels is 40! 
	for (j in 1 : length(block_names)) {
		block <- block_names[[j]]
		print(paste("	block: ", block, sep=""))
		if (!is.null(map[[patient]][[block]])) {
			block_images <- map[[patient]][[block]]
			num_annotations_block <- matrix(data=0, nrow=1, ncol=7)
			for (k in 1 : length(block_images)) {
				print(paste("		image: ", block_images[k], sep=""))
				#each row stores statistics at level of 45*l by 45*l grid, l ranging from 1 to 40, for each annotation/statistic type (40 total also)
				mat_rel_nuclei_fraction <- matrix(data=0, nrow=40, ncol=40)
				fileContents <- read.csv(paste(directory, block_images[k], sep=""), sep=",", header=FALSE) #file is sorted by x and y coordinate already.
				width <- max(fileContents[,1]) + 1 
				height <- max(fileContents[,2]) + 1
				sectionName <- strsplit(block_images[k], ".csv")[[1]][1]
				mat_annot_count <- matrix(data=0, nrow=40, ncol=8)
				hist_all_size_annot <- array(0, dim=c(40,10,8))
				
				for (l in 1 : 40) {
					print(paste("			l: ", l, sep=""))
					#relNucleiFractions_X stores rel nuc fraction averaged across 45 by 45 tiles lying in each 45*l by 45*l tile.  Value of 0 can either mean 45*l by 45*l tile has nonzero number of X annotated 45 by 45 tile, whose rel nuc fraction averages to 0, or, 45*l by 45*l tile contains no X annotated 45 by 45 tiles.  The former case is quite rare (nearly impossible to find)--save for in very late stage necrosis.
					relNucleiFractions_CT <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_NE <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_HBV <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_MVP <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_LE <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_IT <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_PAN <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					relNucleiFractions_annot <- matrix(data=0, nrow=1, ncol=floor(width/(45*l))*floor(height/(45*l)))
					#traverse x and y coordinates (start from (0,0) in the file) in 45*l by 45*l increments for both
 					for (x_coord in seq(0,width,45*l)){ 
						if ((x_coord %% 4500) == 0) {
							print(paste("				x_coord: ", x_coord, sep=""))
						}
						for (y_coord in seq(0,height,45*l)){ 
							cols <- floor(x_coord/45):(floor(x_coord/45)+l-1) + 1 #l cols inclusive.. add 1 because R 1 indexes.
							ind <- cols #stores how to index into rows of file
							p <- 1
							while (p <= l-1) {
								ind <- c(ind, cols + ceiling(width/45)*p)
								p <- p + 1
							}
							#displace each element in ind by appropriate number of rows
							ind <- ind + (floor((y_coord - 0)/45))*ceiling(width/45) #floor((y_coord - 0)/45)) is the "row index"				
							if (length(ind) != 1) {
								print("hmm")
							}
							tiles <- fileContents[ind, ] #issue is here.
							CT <- which(tiles[,3] == "CT")
							NE <- which(tiles[,3] == "NE")
							HBV <- which(tiles[,3] == "HBV")
							MVP <- which(tiles[,3] == "MVP")
							LE <- which(tiles[,3] == "LE")
							IT <- which(tiles[,3] == "IT")
							PAN <- which(tiles[,3] == "PAN")
							annot <- which(tiles[,3]=="CT" || tiles[,3]=="NE" || tiles[,3]=="HBV" || tiles[,3]=="MVP" || tiles[,3]=="LE" || tiles[,3]=="IT" || tiles[,3]=="PAN")
							mat_annot_count[l,] <- c(length(CT), length(NE), length(HBV), length(MVP), length(LE), length(IT), length(PAN), length(annot))
							if (length(CT) == 0) {
								relNucleiFractions_CT[x_coord + (y_coord-1)*floor(width/(45*l))] <-  0
							} else {
								relNucleiFractions_CT[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[CT, 6])
							}
							if (length(NE) == 0) {
								relNucleiFractions_NE[x_coord + (y_coord-1)*floor(width/(45*l))] <-  0
							} else {
								relNucleiFractions_NE[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[NE, 6])
							}
							if (length(HBV) == 0) {
								relNucleiFractions_HBV[x_coord + (y_coord-1)*floor(width/(45*l))] <-  0
							} else {
								relNucleiFractions_HBV[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[HBV, 6])
							}
							if (length(MVP) == 0) {
								relNucleiFractions_MVP[x_coord + (y_coord-1)*floor(width/(45*l))] <- 0
							} else {
								relNucleiFractions_MVP[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[MVP, 6])
							}
							if (length(LE) == 0) {
								relNucleiFractions_LE[x_coord + (y_coord-1)*floor(width/(45*l))] <- 0
							} else {
								relNucleiFractions_LE[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[LE, 6])
							}
							if (length(IT) == 0) {
								relNucleiFractions_IT[x_coord + (y_coord-1)*floor(width/(45*l))] <- 0
							} else {
								relNucleiFractions_IT[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[IT, 6])
							}
							if (length(PAN) == 0) {
								relNucleiFractions_PAN[x_coord + (y_coord-1)*floor(width/(45*l))] <- 0
							} else {
								relNucleiFractions_PAN[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[PAN, 6])
							}
							if (length(annot) == 0) {
								relNucleiFractions_annot[x_coord + (y_coord-1)*floor(width/(45*l))] <- 0
							} else {
								relNucleiFractions_annot[x_coord + (y_coord-1)*floor(width/(45*l))] <- mean(tiles[annot, 6])
							}							
						}
					}

					#counts_X[i] reports number of 45*l by 45*l tiles with rel nuc fraction, averaged across X annotated 45 by 45 tiles, between (i-1)/10 and i/10					
					counts_CT <- matrix(data=0, nrow=1, ncol=10)
					counts_NE <- matrix(data=0, nrow=1, ncol=10)
					counts_HBV <- matrix(data=0, nrow=1, ncol=10)
					counts_MVP <- matrix(data=0, nrow=1, ncol=10)
					counts_LE <- matrix(data=0, nrow=1, ncol=10)
					counts_IT <- matrix(data=0, nrow=1, ncol=10)
					counts_PAN <- matrix(data=0, nrow=1, ncol=10)
					counts_annot <- matrix(data=0, nrow=1, ncol=10)				
				
					for (perc in seq(0.1,1,0.1)) {
						counts_CT[perc*10] <- length(which(relNucleiFractions_CT < perc & relNucleiFractions_CT >= (perc - 0.1)))
						counts_NE[perc*10] <- length(which(relNucleiFractions_NE < perc & relNucleiFractions_NE >= (perc - 0.1)))
						counts_HBV[perc*10] <- length(which(relNucleiFractions_HBV < perc & relNucleiFractions_HBV >= (perc - 0.1)))
						counts_MVP[perc*10] <- length(which(relNucleiFractions_MVP < perc & relNucleiFractions_MVP >= (perc - 0.1)))
						counts_LE[perc*10] <- length(which(relNucleiFractions_LE < perc & relNucleiFractions_LE >= (perc - 0.1)))
						counts_IT[perc*10] <- length(which(relNucleiFractions_IT < perc & relNucleiFractions_IT >= (perc - 0.1)))
						counts_PAN[perc*10] <- length(which(relNucleiFractions_PAN < perc & relNucleiFractions_PAN >= (perc - 0.1)))
						counts_annot[perc*10] <- length(which(relNucleiFractions_annot < perc & relNucleiFractions_annot >= (perc - 0.1)))						
					}
					hist_all_size_annot[l,,1] <- counts_CT
					hist_all_size_annot[l,,2] <- counts_NE
					hist_all_size_annot[l,,3] <- counts_HBV
					hist_all_size_annot[l,,4] <- counts_MVP
					hist_all_size_annot[l,,5] <- counts_LE
					hist_all_size_annot[l,,6] <- counts_IT
					hist_all_size_annot[l,,7] <- counts_PAN
					hist_all_size_annot[l,,8] <- counts_annot		
					
					mat_rel_nuclei_fraction[l, ] <- c(mean(relNucleiFractions_CT), median(relNucleiFractions_CT), median(relNucleiFractions_CT[which(relNucleiFractions_CT <= median(relNucleiFractions_CT))]), median(relNucleiFractions_CT[which(relNucleiFractions_CT > median(relNucleiFractions_CT))]), sd(as.vector(relNucleiFractions_CT)), 
					mean(relNucleiFractions_NE), median(relNucleiFractions_NE), median(relNucleiFractions_NE[which(relNucleiFractions_NE <= median(relNucleiFractions_NE))]), median(relNucleiFractions_NE[which(relNucleiFractions_NE > median(relNucleiFractions_NE))]), sd(as.vector(as.vector(relNucleiFractions_NE))), 
					mean(relNucleiFractions_HBV), median(relNucleiFractions_HBV), median(relNucleiFractions_HBV[which(relNucleiFractions_HBV <= median(relNucleiFractions_HBV))]), median(relNucleiFractions_HBV[which(relNucleiFractions_HBV > median(relNucleiFractions_HBV))]), sd(as.vector(relNucleiFractions_HBV)), 
					mean(relNucleiFractions_MVP), median(relNucleiFractions_MVP), median(relNucleiFractions_MVP[which(relNucleiFractions_MVP <= median(relNucleiFractions_MVP))]), median(relNucleiFractions_MVP[which(relNucleiFractions_MVP > median(relNucleiFractions_MVP))]), sd(as.vector(relNucleiFractions_MVP)), 
					mean(relNucleiFractions_LE), median(relNucleiFractions_LE), median(relNucleiFractions_LE[which(relNucleiFractions_LE <= median(relNucleiFractions_LE))]), median(relNucleiFractions_LE[which(relNucleiFractions_LE > median(relNucleiFractions_LE))]), sd(as.vector(relNucleiFractions_LE)), 
					mean(relNucleiFractions_IT), median(relNucleiFractions_IT), median(relNucleiFractions_IT[which(relNucleiFractions_IT <= median(relNucleiFractions_IT))]), median(relNucleiFractions_IT[which(relNucleiFractions_IT > median(relNucleiFractions_IT))]), sd(as.vector(relNucleiFractions_IT)), 
					mean(relNucleiFractions_PAN), median(relNucleiFractions_PAN), median(relNucleiFractions_PAN[which(relNucleiFractions_PAN <= median(relNucleiFractions_PAN))]), median(relNucleiFractions_PAN[which(relNucleiFractions_PAN > median(relNucleiFractions_PAN))]), sd(as.vector(relNucleiFractions_PAN)), 
					mean(relNucleiFractions), median(relNucleiFractions), median(relNucleiFractions[which(relNucleiFractions <= median(relNucleiFractions))]), median(relNucleiFractions[which(relNucleiFractions > median(relNucleiFractions))]), sd(as.vector(relNucleiFractions)))
				}
				#plot stuff
				colors <- topo.colors(40)
				#for (l in 1:40) {
				for (l in seq(1,40,5)) {				
					if (l == 1) {
						#normalize by number of annotated tiles of size 45*l by 45*l, containing at least one 45*45 CT annotated tile
						plot(seq(0.05, 1, 0.1), hist_all_size_annot[l,, 1]/sum(hist_all_size_annot[l,, 1]), ylim=c(0,1), col=colors[l], type="l", lwd=2)
					} else {
						points(seq(0.05, 1, 0.1), hist_all_size_annot[l,, 1]/sum(hist_all_size_annot[l,, 1]), ylim=c(0,1), col=colors[l], type="l", lwd=2)
					}
				}
				#legend('topright', legend=1:40 , lty=1, col=colors[1:40], bty='n', cex=1.2)
				legend('topright', legend=seq(1,40,5) , lty=1, col=colors[seq(1,40,5)], bty='n', cex=1.2)				
				
				mat_rel_nuclei_fraction[which(is.na(mat_rel_nuclei_fraction))] <- 0
				colnames(mat_rel_nuclei_fraction) <- c("CT_mean", "CT_median", "CT_Q1", "CT_Q3", "CT_sd", "NE_mean", "NE_median", "NE_Q1", "NE_Q3", "NE_sd", "HBV_mean", "HBV_median", "HBV_Q1", "HBV_Q3", "HBV_sd", "MVP_mean", "MVP_median", "MVP_Q1", "MVP_Q3", "MVP_sd", "LE_mean", "LE_median", "LE_Q1", "LE_Q3", "LE_sd", "IT_mean", "IT_median", "IT_Q1", "IT_Q3", "IT_sd", "PAN_mean", "PAN_median", "PAN_Q1", "PAN_Q3", "PAN_sd", "all_annot_mean", "all_annot_median", "all_annot_Q1", "all_annot_Q3", "all_annot_sd")
				write.table(mat_rel_nuclei_fraction, paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/RelNucFraction/", sectionName, ".csv", sep=""), col.names=TRUE, row.names=FALSE)				
				
				#Annotation percentages
				len_CT <- length(which(fileContents[,3] == "CT"))
				len_NE <- length(which(fileContents[,3] ==  "NE"))
				len_HBV <- length(which(fileContents[,3] ==  "HBV"))
				len_MVP <- length(which(fileContents[,3] ==  "MVP"))
				len_LE <- length(which(fileContents[,3] ==  "LE"))
				len_IT <- length(which(fileContents[,3] ==  "IT"))
				len_PAN <- length(which(fileContents[,3] ==  "PAN"))
				
				vec <- t(c(len_CT, len_NE, len_HBV, len_MVP, len_LE, len_IT, len_PAN))
				colnames(vec) <- c("CT", "NE", "HBV", "MVP", "LE", "IT", "PAN")
				write.table(vec, paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""), append=TRUE, row.names=FALSE, col.names=TRUE)
				write.table(vec/sum(vec), paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""), append=TRUE, row.names=FALSE, col.names=FALSE)
				
				num_annotations_block[1,1] <- num_annotations_block[1,1] + len_CT
				num_annotations_block[1,2] <- num_annotations_block[1,2] + len_NE
				num_annotations_block[1,3] <- num_annotations_block[1,3] + len_HBV
				num_annotations_block[1,4] <- num_annotations_block[1,4] + len_MVP
				num_annotations_block[1,5] <- num_annotations_block[1,5] + len_LE
				num_annotations_block[1,6] <- num_annotations_block[1,6] + len_IT
				num_annotations_block[1,7] <- num_annotations_block[1,7] + len_PAN

				num_annotations_patient[1,1] <- num_annotations_patient[1,1] + len_CT
				num_annotations_patient[1,2] <- num_annotations_patient[1,2] + len_NE
				num_annotations_patient[1,3] <- num_annotations_patient[1,3] + len_HBV
				num_annotations_patient[1,4] <- num_annotations_patient[1,4] + len_MVP
				num_annotations_patient[1,5] <- num_annotations_patient[1,5] + len_LE
				num_annotations_patient[1,6] <- num_annotations_patient[1,6] + len_IT
				num_annotations_patient[1,7] <- num_annotations_patient[1,7] + len_PAN

				z_level <- as.numeric(substr(strsplit(sectionName, "_p")[[1]][2], 1, 2))
				num_annotations_z[z_level,1] <- num_annotations_z[z_level,1] + len_CT
				num_annotations_z[z_level,2] <- num_annotations_z[z_level,2] + len_NE
				num_annotations_z[z_level,3] <- num_annotations_z[z_level,3] + len_HBV
				num_annotations_z[z_level,4] <- num_annotations_z[z_level,4] + len_MVP
				num_annotations_z[z_level,5] <- num_annotations_z[z_level,5] + len_LE
				num_annotations_z[z_level,6] <- num_annotations_z[z_level,6] + len_IT
				num_annotations_z[z_level,7] <- num_annotations_z[z_level,7] + len_PAN
			}
			for (k in 1 : length(block_images)) {
				sectionName <- strsplit(block_images[k], ".csv")[[1]][1]
				vec <- read.table(paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""))
				z_level <- as.numeric(substr(strsplit(sectionName, "_p")[[1]][2], 1, 2))
				write.table(vec[1,]/sum(num_annotations_z[z_level,]), paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""), append=TRUE, row.names=FALSE, col.names=FALSE)
				write.table(vec[1,]/sum(num_annotations_patient[1,]), paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""), append=TRUE, row.names=FALSE, col.names=FALSE)
			}
		}
	}
	images <- map2[[patient]]
	for (img in 1 : length(images)) {
		sectionName <- strsplit(images[img], ".csv")[[1]][1]
		vec <- read.table(paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""))
		write.table(vec[1,]/sum(num_annotations_block[1,]), paste("C:/Users/Nameeta/Chantal/results/RelNucFraction_Annotations/", patient, "/PercAnnotations/", sectionName, ".csv", sep=""), append=TRUE, row.names=FALSE, col.names=FALSE)
	}	
}
