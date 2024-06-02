### compute structural score (spatial_chaos) for each ion image of an MSI 

library(mmand)
## dilation_and_erosion: apply morphological operations to a binary image

# im: 2D matrix of a binary image

dilation_and_erosion <- function(im){
  dilate_mask <- rbind(c(0,1,0),c(1,1,1),c(0,1,0))
  erode_mask <- rbind(c(1,1,1),c(1,1,1),c(1,1,1))
  
  return(erode(dilate(im, dilate_mask), erode_mask))
}

## count_func: count the number of separate objects in a binary image

# im: 2D matrix of a binary image

count_func <- function(im){
  k <- shapeKernel(c(3,3), type="diamond") # 4 connectivity
  Clumps <- components(im, k)
  
  if(length(which(Clumps > 0)) == 0){ # otherwise it will return -inf 
    return(0)
  }
  else{
    tot <- max(Clumps, na.rm=TRUE)
    return(tot)
  }
}

## counts_levels: Created binary images

# im: The values of intensities are between 0 and 1
# nlevels: integer specifying the number of binary image to create from the original

count_levels <- function(im,nlevels){
  
  # compute for each level where to cut
  levels <- seq(from=0,to=1 ,by=1/nlevels)[-(nlevels+1)]
  
  nb_tot_objects <- 0
  nb_tot_px <- 0
  
  for(i in 1:length(levels)){
    # create the binary image 
    bw <- (im > levels[i])
    #bw <- dilation_and_erosion(bw)
    #image(bw,col = c("black", "white"))
    nb_px<- length(which(bw >0))
    
    # comute the number of object in each binary image
    if(nb_px > 1){
        object_num = sum(count_func(bw))
        
        if(nb_px <= object_num){
            nb_tot_objects <- nb_tot_objects +1 
            nb_tot_px <-  nb_tot_px  +1
        }
        else{
            nb_tot_objects <- nb_tot_objects + object_num -2
            nb_tot_px <-  nb_tot_px +  nb_px
        }
            
    }
    if(nb_px <= 1){
      nb_tot_objects <- nb_tot_objects +1 
      nb_tot_px <-  nb_tot_px  +1
    }
    
  }
  # return a value between 0 and 1 for each ionic image
  m_chaos <- 1-  (nb_tot_objects/nb_tot_px)
  return(m_chaos)
}

## measure_spatial_chaos: compute the value of chaos from an ionic image

# im: The values of intensities are between 0 and 1
# nlevels: integer specifying the number of binary image to create from the original


measure_spatial_chaos <- function(im,nlevels){
  if( sum(im,na.rm = TRUE) <= 0 ){return(NA)}
  im_clean <- im/max(im,na.rm=TRUE)
  chaos_val  <- count_levels( im_clean,nlevels)
  
  return(chaos_val)
}

hot_spot_removal <- function(im,q){
  xic_q <- quantile(im,q)
  im[im > xic_q] <- xic_q 
  return(im)
}

make_mat_image <- function(im,crd){
  mat <- matrix(data=0,nrow=max(crd[,2]),ncol=max(crd[,1]))
  
  for(i in 1:length(crd[,1])){
    mat[crd[i,2],crd[i,1]] <- im[i]
  }
  return(mat)
}

tic_norm <- function(msi){
  msi[is.na(msi)] <- 0
  msi_norm <- msi
  for(i in 1:ncol(msi)){
    msi_norm[,i] <- msi[,i]/sum(msi[,i])
  }
  return(msi_norm)
}

compute_structure <- function(folder_list,bin){
  for(f in 1:length(folder_list)){
    outputfolder <- folder_list[f]
    crd <- read.csv(paste(outputfolder,"crd.csv",sep=""), header=FALSE)
    crd[,1] <- crd[,1]+1
    crd[,2] <- crd[,2]+1
    msi <- read.csv(paste(outputfolder,"msi.csv",sep=""), header=FALSE)
    msi <- as.matrix(msi)
    #msi <- tic_norm(msi)
    chaos_vec <- rep(0,nrow(msi))
    for(i in 1:nrow(msi)){
      imd <- msi[i,]
      #imd <- hot_spot_removal(imd,0.99)
      mat <- make_mat_image(imd,crd)
      chaos_vec[i] <- measure_spatial_chaos(mat,bin)
      print(i)
    }
    write.table(chaos_vec,file=paste(outputfolder,"_spatial_score.csv",sep=""),col.names=F,row.names=F)
    #return(chaos_vec)
  }
}


# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("One argument must be supplied.")
}

# Access the arguments
arg1 <- args[1]


compute_structure(paste(arg1,"_",sep=""),30)