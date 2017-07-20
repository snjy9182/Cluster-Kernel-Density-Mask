trackll <- createTrackll(interact = T, cores = 5)

trackl <- trackll[[3]]
trackl.m <- mergeTrack(trackl)

dens <- calculateKernelDensity(trackl.m)
cluster.region <- getClusters(trackl.m, dens, p = 0.3)

trackl.masked <- clusterMask(trackl, cluster.region)
trackl.m.masked <- mergeTrack(trackl.masked)
plotTrackl(trackl.masked)


# Merge into one list
################################################################################################
library(dplyr)
mergeTrack = function(track.list, plot = T){
  if (length(track.list[[1]]) == 3){
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z")]
  } else {
    df <- bind_rows(track.list, .id = "Trajectory")[, c("x", "y", "z", "Frame")]
  }
  if (plot){
    plot(df[[1]], df[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", cex = .2)
  }
  return (df);
}

# Extract contours
################################################################################################

library(MASS)
library(sp)

calculateKernelDensity = function(trackl.merged){
  dens <- MASS::kde2d(trackl.merged[[1]], trackl.merged[[2]], n=200, lims=c(c(0, 128), c(0, 128)))
  return(dens);
}

getClusters = function (trackl.merged, dens, p = 0.3, num.clusters = -1){
  x = trackl.merged[[1]]
  y = trackl.merged[[2]]
  # the contours to plot
  prob <- c(p)
  dx <- diff(dens$x[1:2])
  dy <- diff(dens$y[1:2])
  sz <- sort(dens$z)
  c1 <- cumsum(sz) * dx * dy 
  levels <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
  })
  
  ## Points within polygons
  ls <- contourLines(dens, level=levels)
  
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
  nuclei <- list()
  for (i in 1:length(ls)){
    nuclei[[i]] <- point.in.polygon(x, y, ls[[i]]$x, ls[[i]]$y)
  }
  ## Plot
  trackl.merged$region <- factor(Reduce("+", nuclei))
  plot(y ~ x, col=region, data=trackl.merged, cex = 0.2)
  contour(dens, levels=levels, labels=prob, add=T)
  return (trackl.merged$region)
}

# Mask trackl
################################################################################################

clusterMask = function (track.list, cluster.region){
  masked.track.list = list();
  index.mask = 1;
  index = 1;
  for(i in 1:length(track.list)){
    mask = TRUE;
    for (j in 1:nrow(track.list[[i]])){
      if (cluster.region[[index]] == 1){
        mask = FALSE;
      }
      index = index + 1;
    }
    if (!mask){
      masked.track.list[index.mask] <- track.list[i];
      index.mask = index.mask + 1;
    }
  }
  return (masked.track.list);
}

# Plotting trackl
################################################################################################


plotTrackLines = function(track.list){
  plot(track.list[[3]][[1]], track.list[[3]][[2]], type = "l", xlim = c(0, 128), ylim = c(0, 128))
  for(i in 4:length(track.list)){
    lines(track.list[[i]][[1]], track.list[[i]][[2]])
  }
}

