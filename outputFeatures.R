library(raster)
require('igraph')
path <- "C:/Users/Nameeta/Chantal/data/W12-1-1-D.2"
annotation_files <- list.files(paste(path, "/annotations", sep=""))
scalev <- 45

for (file_ind in 1:length(annotation_files)) {
	print(file_ind)
	#do MST stuff
	annotation_file <- annotation_files[file_ind] 
	sliceName <- strsplit(annotation_file, ".xml")
	if (nchar(sliceName) < nchar(annotation_file)) {
		nucleiStatistics <- read.csv(paste(path, "/NucleiStatistics/", sliceName, ".v2.csv", sep=""), sep=";")
		width <- max(nucleiStatistics[,2]) + 45
		height <- max(nucleiStatistics[,3]) + 45
		num_width <- ceiling(width/scalev)
		num_height <- ceiling(height/scalev)
		output_matrix <- matrix(data=0, nrow=num_width*num_height, ncol=14)
		
		#Brightness/rel frac area features 	
		for (tile_ind in 1:nrow(nucleiStatistics)) {
			innerX <- nucleiStatistics[tile_ind,2]
			innerY <- nucleiStatistics[tile_ind,3]
			if (innerX == 0) {
				innerX = 1;
			}
			if (innerY == 0) {
				innerY = 1;
			}
			relAreaNuclei <- nucleiStatistics[tile_ind, 9]	
			nucleiBrightness <-  nucleiStatistics[tile_ind,10]*relAreaNuclei
			backgroundBrightness <- nucleiStatistics[tile_ind, 11]*(1-relAreaNuclei)
			output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),10] <- relAreaNuclei		
			output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),11] <- nucleiBrightness
			output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),12] <- backgroundBrightness
			output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),13] <- output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),11] + output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),12]
			#add ratio of brightness
			if (backgroundBrightness > 0 && nucleiBrightness > 0) {
				output_matrix[ceiling(innerY/scalev)*num_width + ceiling(innerX/scalev),14] <- log(nucleiBrightness / backgroundBrightness)
			}
		}
		
		nucleiCoordinates <- read.csv(paste(path, "/NucleiCoordinates/", sliceName, ".v2.csv", sep=""), sep=";") 		 
		coordinates <- as.matrix(nucleiCoordinates[,c(2,3)])	
		num_nuclei <- nrow(coordinates)		
		names(coordinates) <- c("X", "Y")
		#map each nuclei to a tile
		max_capacity <- round(500*(num_nuclei/(num_width*num_height))) #tweak this value based on distribution of nuclei across the regions.		
		nuclei_all_region = array(0, dim=c(num_height,num_width,max_capacity,2))	
		num_nuclei_all_region <- matrix(data=0, nrow=num_height, ncol=num_width)
		for (k in 1:num_nuclei) {
			location <- coordinates[k,]
			row_ind <- ceiling((as.numeric(location[2])/height)*num_height) #as.numeric necessary?
			col_ind <- ceiling((as.numeric(location[1])/width)*num_width)
			num_nuclei_region <- num_nuclei_all_region[row_ind, col_ind]
			nuclei_all_region[row_ind, col_ind, num_nuclei_region + 1,1] <- location[1]
			nuclei_all_region[row_ind, col_ind, num_nuclei_region + 1,2] <- location[2]
			num_nuclei_all_region[row_ind, col_ind] <- num_nuclei_region + 1
		}
		
		#MST features
		#MSTs_tiles <- as.list(numeric(num_height*num_width))
		#dim(MSTs_tiles) <- c(num_height, num_width)				
		for (i in 1:num_width) {
			for (j in 1:num_height) {
				num_nuclei_region <- num_nuclei_all_region[j, i]
				nuclei <- (nuclei_all_region[j,i,,])[1:num_nuclei_region,]
				if (num_nuclei_region > 1) { # build matrix representing fully connected graph with no self loops
					mat <- matrix(1, num_nuclei_region, num_nuclei_region)
					g <- graph.adjacency(mat) # build graph object
					E(g)$weight <- as.matrix(dist(nuclei)) # lists edge weights in row order
					
					#find angle of each edge with respect to positive x unit vector
					angle_mat <- matrix(-1,nrow=num_nuclei_region, ncol=num_nuclei_region)
					diff_x <- outer(nuclei[,1],nuclei[,1],"-")
					diff_y <- outer(nuclei[,2],nuclei[,2],"-")	
					
					ind1 <- which(diff_x == 0 & diff_y < 0)
					ind2 <- which(diff_x == 0 & diff_y >= 0)
					ind3 <- which(diff_x > 0 & diff_y >= 0)
					ind4 <- which(diff_x < 0 & diff_y > 0)
					ind5 <- which(diff_x < 0 & diff_y <= 0)	
					ind6 <- which(diff_x > 0 & diff_y < 0)					
					
					angle_mat[ind1] <- (3*pi)/2 
					angle_mat[ind2] <- (pi)/2 			
					angle_mat[ind3] <- atan(diff_y[ind3]/diff_x[ind3])
					angle_mat[ind4] <- pi - abs(atan(diff_y[ind4]/diff_x[ind4]))
					angle_mat[ind5] <- atan(diff_y[ind5]/diff_x[ind5]) + 180
					angle_mat[ind6] <- 2*pi - atan(diff_y[ind6]/diff_x[ind6])

					stopifnot(length(which(angle_mat < 0)) == 0)
					E(g)$angle <- angle_mat
					output_matrix[(j-1)*num_width + i,6] <- num_nuclei_region
					if (num_nuclei_region > 1) {
						output_matrix[(j-1)*num_width + i,7] <- mean(E(g)$weight)
						output_matrix[(j-1)*num_width + i,8] <-	max(E(g)$weight)				
						output_matrix[(j-1)*num_width + i,9] <-	mean(E(g)$angle)
					}					
				} 
			}
		}
		
		#Neighbor features and annotation
		#pass width, height, annotation name to perl
		setwd('C:/strawberry/perl/bin')
		system(paste("perl C:/Users/Nameeta/Chantal/src/Perl/processAnnotations.pl ", width, height, sliceName, 1)) #just changed for now
		setwd('C:/Users/Nameeta/Chantal/src/R')
		#read in annotations file and make box plot of comparison of features.		
		annotationsOverlay <- read.csv(paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Annotations/", sliceName, ".csv", sep=""), sep=",")
		CT_tiles <- annotationsOverlay[which(annotationsOverlay[,3]=="CT"),1:2]
		NE_tiles <- annotationsOverlay[which(annotationsOverlay[,3]=="NE"),1:2]
		PAN_tiles <- annotationsOverlay[which(annotationsOverlay[,3]=="PAN"),1:2]
		HBV_tiles <- annotationsOverlay[which(annotationsOverlay[,3]=="HBV"),1:2]		
		
		for (i in 1:num_width) {
			for (j in 1:num_height) {		
				if (i > 5 && j > 5 && (i < (num_width - 5)) && (j < (num_height - 5))) {
					nuclei_area <- output_matrix[(j-1)*num_width + i, 10]	
					nuclei_brightness <- output_matrix[(j-1)*num_width + i, 11]
					background_brightness <- output_matrix[(j-1)*num_width + i, 12]
					total_brightness <- output_matrix[(j-1)*num_width + i, 13]						
					ind_vec <- c()
					for (a in -5:5) {
						for (b in -5:5) {
							ind_vec <- c(ind_vec, (j-a)*num_width + i-b)
						}
					}
					mean_nuclei_area <- mean(output_matrix[ind_vec, 10])
					mean_nuclei_brightness <- mean(output_matrix[ind_vec, 11])
					mean_background_brightness <- mean(output_matrix[ind_vec, 12])
					mean_total_brightness <- mean(output_matrix[ind_vec, 13])
					if (nuclei_area != 0 || mean_nuclei_area != 0) {
						output_matrix[(j-1)*num_width + i,2] <-  min(nuclei_area, mean_nuclei_area)/max(nuclei_area, mean_nuclei_area)
					}
					if (nuclei_brightness != 0 || mean_nuclei_brightness != 0) {
						output_matrix[(j-1)*num_width + i,3] <- min(nuclei_brightness, mean_nuclei_brightness)/max(nuclei_brightness, mean_nuclei_brightness)
					}
					if (background_brightness != 0 || mean_background_brightness != 0) {
						output_matrix[(j-1)*num_width + i,4] <- min(background_brightness, mean_background_brightness)/max(background_brightness, mean_background_brightness)
					}
					if (total_brightness != 0 || mean_total_brightness != 0) {
						output_matrix[(j-1)*num_width + i,5] <- min(total_brightness, mean_total_brightness)/max(total_brightness, mean_total_brightness)
					}
				}
				if (output_matrix[(j-1)*num_width + i,1] > 0) {
					print(paste("WTF: ", output_matrix[(j-1)*num_width + i,1]), sep="")
				}
			
				if (length(intersect(which(CT_tiles[,1]==i), which(CT_tiles[,2]==j)) > 0)) {
					output_matrix[(j-1)*num_width + i,1] <- 1
				} else if (length(intersect(which(NE_tiles[,1]==i), which(NE_tiles[,2]==j)) > 0)) {
					output_matrix[(j-1)*num_width + i,1] <- 2					
				} else if (length(intersect(which(PAN_tiles[,1]==i), which(PAN_tiles[,2]==j)) > 0)) {
					output_matrix[(j-1)*num_width + i,1] <- 3					
				} else if (length(intersect(which(HBV_tiles[,1]==i), which(HBV_tiles[,2]==j)) > 0)) {
					output_matrix[(j-1)*num_width + i,1] <- 4					
				} else {
					output_matrix[(j-1)*num_width + i,1] <- 5	
				}
			}		
		}
		
		normalized_output_matrix <- output_matrix
		for (j in 2:ncol(normalized_output_matrix)) {
			mn <- mean(normalized_output_matrix[,j])
			sdv <- sd(normalized_output_matrix[,j])
			normalized_output_matrix[,j] <- (normalized_output_matrix[,j] - mn)/sdv
			if((mean(normalized_output_matrix[,j]) < 0.000001) && (abs(sd(normalized_output_matrix[,j])-1) < 0.000001))  {
				print (paste("sensible: ", j))
			}
		}
		
		#annotated matrix
		annotated_output_matrix <- output_matrix[which(output_matrix[,1] < 5),]
		normalized_annotated_output_matrix <- annotated_output_matrix
		for (j in 2:ncol(normalized_annotated_output_matrix)) {
			mn <- mean(normalized_annotated_output_matrix[,j])
			sdv <- sd(normalized_annotated_output_matrix[,j])
			normalized_annotated_output_matrix[,j] <- (normalized_annotated_output_matrix[,j] - mn)/sdv
			if((mean(normalized_annotated_output_matrix[,j]) < 0.000001) && (abs(sd(normalized_annotated_output_matrix[,j])-1) < 0.000001))  {
				print (paste("sensible: ", j))
			}
		}		
		
		# #choose rows of annotated matrix so ROIs equally represented
		# CT_inds <- which(output_matrix[,1]==1)
		# HBV_inds <- which(output_matrix[,1]==4)
		# whc <- which(output_matrix[,1]==2)
		# NE_inds <- whc[sample(1:length(whc), 500, replace=F)]
		# equalized_annotated_output_matrix <- output_matrix[c(CT_inds, HBV_inds, NE_inds),]
		
		# normalized_equalized_annotated_output_matrix <- equalized_annotated_output_matrix
		# for (j in 2:ncol(normalized_equalized_annotated_output_matrix)) {
			# mn <- mean(normalized_equalized_annotated_output_matrix[,j])
			# sdv <- sd(normalized_equalized_annotated_output_matrix[,j])
			# normalized_equalized_annotated_output_matrix[,j] <- (normalized_equalized_annotated_output_matrix[,j] - mn)/sdv
			# if((mean(normalized_equalized_annotated_output_matrix[,j]) < 0.000001) && (abs(sd(normalized_equalized_annotated_output_matrix[,j])-1) < 0.000001))  {
				# print (paste("sensible: ", j))
			# }
		# }		
		
		write.table(output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_all_", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE) 
		write.table(normalized_output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_normalized_all_", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE)
		write.table(annotated_output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_annotated_", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE)
		write.table(normalized_annotated_output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_normalized_annotated_", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE)		
		#write.table(equalized_annotated_output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_annotated_equalized", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE)
		#write.table(normalized_equalized_annotated_output_matrix, paste("C:/Users/Nameeta/Chantal/results/W12-1-1-D.2/Features/feat_normalized_annotated_equalized", sliceName, ".txt", sep=""), row.names=FALSE, col.names=FALSE)		
	}
}	

