tl <- trackll.linked[[10]]
num.cells = 10
# Merge into one list
################################################################################################
df <- NULL
temp <- tl
for (i in 1:length(temp)){
  temp[[i]] <- subset(temp[[i]], select = -c(3, 4))
  df <- rbind(df, temp[[i]])
}
#plot(df[[1]], df[[2]], xlim = c(0, 128), ylim = c(0, 128), xlab = "x", ylab = "y", cex = .2)
# Extract contours
################################################################################################

library(MASS)
x = df[[1]]
y = df[[2]]
dens <- MASS::kde2d(x, y, n=200, lims=c(c(0, 128), c(0, 128)))  # don't clip the contour

# the contours to plot
prob <- c(0.6)
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
while (length(ls) > num.cells){
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
  
nuclei <- list()
for (i in 1:length(ls)){
  nuclei[[i]] <- point.in.polygon(x, y, ls[[i]]$x, ls[[i]]$y)
}
## Plot
df$region <- factor(Reduce("+", nuclei))
plot(y ~ x, col=region, data=df, cex = 0.2)
contour(dens, levels=levels, labels=prob, add=T)
