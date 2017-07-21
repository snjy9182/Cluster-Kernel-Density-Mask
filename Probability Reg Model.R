file = file.choose()
data = read.csv(file)

average.track.length = data$Points/data$Trajectories
probability = data$Probabillity
summary(fit <- lm(probability ~ average.track.length))


plot(probability ~ average.track.length, xlab = "Average Track Length (Unfiltered)", ylab = "Kernel Density Probability", main = "Kernel Density Probability Estimation from Variable Track Lists")
abline(fit)

