#### maskTrackll.R
#### Wu Lab, Johns Hopkins University
#### Author: Sun Jay Yoo
#### Date: July 21, 2017

## maskTrackll-methods
##
##
###############################################################################
##' @name maskTrackll
##' @aliases maskTrackll
##' @title maskTrackll
##' @rdname maskTrackll-methods
##' @docType methods
##'
##' @description mask track lists and lists of track lists using kernel density clusters

##' @usage 
##' maskTrackll(trackll, p = NULL)
##' 
##' maskTrackl(track.list, mask)
##' 
##' createMask(track.list, kernel.density, p = NULL, num.clusters = -1, plot = T)
##' 
##' kernelDensity(track.list)
##' 
##' plotTrackPoints(track.list)
##' 
##' plotTrackLines(track.list)
##' 

##' @param trackll An uncensored/unfiltered list of track lists.
##' @param track.list A single uncensored/filtered track list.
##' @param p Kernel density probability (calculated automatically by default or given through user input).
##' @param mask A binary mask created using createMask().
##' @param kernel.density A two dimensional kernel density estimation of all track coordinates.
##' @param num.clusters Maximum number of clusters to create (automatic by default).
##' @param plot Plot all track coordinates with the mask (in red) and contour lines.

##' @details
##' The general method for creating a masked track list from a single track list begins by 
##' first calculating its kernel density using kernelDensity(), 
##' creating the mask using this value createMask() (while adjusting the kernel density probability [p] as needed), 
##' then generating the masked track list using maskTrackl. The reason why these three steps are separated is
##' to allow for quick repeated adjustments to the kernel density probability (p), as the other two steps can take more time.
##' 
##' The value for the kernel density probability (p) is automatically calculated by default using 
##' a regression model estimating the approximate desired probability using the track list's average track length.
##' Thus, it is highly recommended to use uncensored/unfiltered track lists.
##' 
##' In order to calculate a masked list of track lists, use maskTrackll(). 
##' All three steps of the masking process are combined and automated across track lists.
##' By default, the kernel density probabilty is automatically calculated, but can be globally overrided by providing user input.
##' 
##' plotTrackPoints() and plotTrackLines() are useful tools to visualize single track lists at any point in the analysis process.

##' @examples
##' ## Process for masking a single track list
##' 
##' #Calculate kernel density
##' kd <- kernelDensity(trackl)
##' 
##' #Create binary mask automatically
##' mask <- createMask(trackl, kd)
##' 
##' #Adjust binary mask manually
##' mask <- createMask(trackl, kd, p = 0.3)
##' mask <- createMask(trackl, kd, p = 0.4)
##' mask <- createMask(trackl, kd, p = 0.45)
##' 
##' #Generate the masked track list
##' masked.trackl <- maskTrackl(trackl, mask)
##' 
##' #Plot track list
##' plotTrackPoints(masked.trackl)
##' plotTrackLines(masked.trackl)
##' 
##' ## Process for making lists of track lists
##' 
##' #Calculate p automatically
##' masked.trackll <- maskTrackll(trackll)
##' 
##' #Use a global p value (not recommended), useful if track lists are uniformly distributed
##' masked.trackll <- maskTrackll(trackll, p = 0.3)

##' @export kernelDensity
##' @export createMask 
##' @export maskTrackl
##' @export maskTrackll
##' @export plotTrackPoints
##' @export plotTrackLines

##' @importFrom dplyr bin_rows
##' @importFrom dplyr MASS::kde2d
##' @importFrom sp point.in.polygon

###############################################################################

#library(dplyr) #bind_rows, MASS::kde2d
#library(sp) #point.in.polygon

#### kernelDensity ####

#Returns kernel density

kernelDensity = function (track.list){
  
    #Merge track list into a single dataframe
    df <- mergeTracks(track.list)
    
    #Calculate kernel density from dataframe
    dens <- MASS::kde2d(df[[1]], df[[2]], n=200, lims=c(c(0, 128), c(0, 128)));

  	return (dens);
}

#### createMask ####

#Returns binary mask and plots

createMask = function (track.list, kernel.density, p = NULL, num.clusters = -1, plot = T){
	#Store all merged track coordinate points into a dataframe
	df <- mergeTracks(track.list)
	
	if (is.null(p)){
	  p = -0.1207484 + 0.3468734*(nrow(df)/length(track.list))
	}
	
	if (p <= 0 || p >= 1){
	  cat("ERROR: Need valid probability (p) or automatic calculation is not valid.")
	}

	# Calculate contours to plot
	prob <- c(p)
	dx <- diff(kernel.density$x[1:2])
	dy <- diff(kernel.density$y[1:2])
	sz <- sort(kernel.density$z)
	c1 <- cumsum(sz) * dx * dy 
	levels <- sapply(prob, function(x){ 
		approx(c1, sz, xout = 1 - x)$y
	})

	#Create the countour polygon with using coordinate points
	ls <- contourLines(kernel.density, level=levels)

	#Keep only the largest user-specified number of clusters, if given
	if (num.clusters > 0){
		while (length(ls) > num.clusters){
		 	noise = 0;
		  	min = Inf;
		  	for (i in 1:length(ls)){
		    	if(length(ls[[i]][[2]]) < min){
		      		noise = i
		      		min = length(ls[[i]][[2]])
		    	}
		  	}
		  	ls[[noise]] <- NULL
		}
	}

	#Use coordinate coordinate polygon to create the cluster shape
	cluster <- list()
	for (i in 1:length(ls)){
		cluster[[i]] <- point.in.polygon(df[[1]], df[[2]], ls[[i]]$x, ls[[i]]$y)
	}

	#Create binary mask of track coordinates
	df$region <- factor(Reduce("+", cluster))

	#Plot with mask and contour
	if(plot){
	  title = paste(getTrackFileName(track.list),"Mask with Kernel Density Probability (p) of", round(p, digits = 3), sep = " ");
		plot(df[[2]] ~ df[[1]], col=region, data=df, xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", main = title, cex = .1)
		contour(kernel.density, levels=levels, labels=prob, add=T)
	}
	cat("\n", getTrackFileName(track.list), "masked at a kernel density probability of", round(p, digits = 3), ".\n")
  return(df$region)
}

#### maskTrackll ####

#Creates masked track list

maskTrackl = function(track.list, mask){
  #Instantiate a masked track list with indexing variables
  masked.track.list = list();
  masked.track.list.names = list();
  index.mask = 1;
  index = 1;
  
  #Loop through all tracks
  for(i in 1:length(track.list)){
    mask.bool = TRUE;
    
    #Remove any tracks outside mask
    for (j in 1:nrow(track.list[[i]])){
      if (mask[[index]] == 1){
        mask.bool = FALSE;
      }
      index = index + 1;
    }
    if (!mask.bool){
      masked.track.list[index.mask] <- track.list[i];
      index.mask = index.mask + 1;
      masked.track.list.names[1 + length(masked.track.list.names)] = names(track.list[i]);
    }		
  }
  names(masked.track.list) <- masked.track.list.names;
  #Return masked track list
  return (masked.track.list);
}

#### mergeTracks ####

mergeTracks = function(track.list){
  if (length(track.list[[1]]) == 3){
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z")]
  } else {
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z", "Frame")]
  }
  
  return (df);
}

#### maskTrackll ####

maskTrackll = function (trackll, p = NULL){
  masked.trackll <- list()
  for (i in 1:length(trackll)){
    kd = kernelDensity(trackll[[i]]);
    mask = createMask(trackll[[i]], kd, p = p)
    masked.trackll[[i]] <- maskTrackl(trackll[[i]], mask)
  }
  names(masked.trackll) <- names(trackll)
  cat("\nAll tracks lists masked.\n")
  return(masked.trackll)
}

#### plotTrackPoints ####

plotTrackPoints = function(track.list){
  df <- mergeTracks(track.list)
  plot(df[[1]], df[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", main = getTrackFileName(track.list), cex = .1);
}

#### plotTrackLines ####

plotTrackLines = function(track.list){
  plot(track.list[[1]][[1]], track.list[[1]][[2]], type = "l", xlim = c(0, 128), ylim = c(0, 128), main = getTrackFileName(track.list))
  for(i in 2:length(track.list)){
    lines(track.list[[i]][[1]], track.list[[i]][[2]])
  }
}

