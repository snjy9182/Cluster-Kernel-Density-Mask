# Collecting Data
################################################################################################
trackl <- .readDiaSessions(interact = T)
trackl.linked <- .linkSkippedFrames(trackl, tolerance = 5, maxSkip = 10)

# Filter
################################################################################################
tl <- list()
name.list <- list()
counter = 1
for (i in 1:length(trackl.linked)){
    if (nrow(trackl.linked[[i]]) >= 2){
        tl[[counter]] <- trackl.linked[[i]]
        name.list[[counter]] <- names(trackl.linked[i])
        counter = counter + 1
    }
}
names(tl) <- name.list

# Merge into one list
################################################################################################
df <- NULL
temp <- tl
for (i in 1:length(temp)){
    temp[[i]] <- subset(temp[[i]], select = -c(3, 4))
    df <- rbind(df, temp[[i]])
}
plot(df[[1]], df[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", cex = .2)

# Merge into average list
################################################################################################
df.avg <- NULL
temp <- tl
for (i in 1:length(temp)){
    temp[[i]] <- data.frame(x = mean(temp[[i]][[1]]), y = mean(temp[[i]][[2]]))
    df.avg <- rbind(df.avg, temp[[i]])
}
plot(df.avg[[1]], df.avg[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", cex = .2)

# Count points within a square
################################################################################################
cl <- dbscan(df, eps = 0.69, minPts = 1000)
hullplot(df, cl = dbscan(df, eps = 6, minPts = 500))
         
count = 0
for (i in 1:length(df[[1]])){
    if (df[[1]][[i]] > 20 && df[[1]][[i]] < 40 && df[[2]][[i]] > 85 && df[[2]][[i]] < 100){
        count = count + 1
    }
}

# Plot contours
################################################################################################

x = df[[1]]
y = df[[2]]
subplot(
    plot_ly(x = x, type = "histogram"),
    plotly_empty(),
    plot_ly(x = x, y = y, type = "histogram2dcontour"),
    plot_ly(y = y, type = "histogram"),
    nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
    shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
)

# Extract contours
################################################################################################

library(MASS)
x = df[[1]]
y = df[[2]]
dens <- MASS::kde2d(x, y, n=200, lims=c(c(0, 128), c(0, 128)))  # don't clip the contour

# the contours to plot
prob <- c(0.4)
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})

## Points within polygons
library(sp)
ls <- contourLines(dens, level=levels)
nuclei3 <- point.in.polygon(x, y, ls[[3]]$x, ls[[3]]$y)
nuclei2 <- point.in.polygon(x, y, ls[[2]]$x, ls[[2]]$y)
nuclei1 <- point.in.polygon(x, y, ls[[1]]$x, ls[[1]]$y)

## Plot
df$region <- factor(nuclei1 + nuclei2 + nuclei3)
plot(y ~ x, col=region, data=df, cex = 0.2)
contour(dens, levels=levels, labels=prob, add=T)
