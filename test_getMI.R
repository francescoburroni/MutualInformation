# Example usage:
set.seed(1310)
x <- rnorm(1000)
y <- rnorm(1000)

source("getMI.R")
MI_nats <- getMI(x, y, nBins=20, smoothingValue = 0.5, units = "nats")
MI_bits <- getMI(x, y, nBins=20, smoothingValue = 0.5, units = "bits")
cat(sprintf("Mutual Information: %.4f nats\n", MI_nats))
cat(sprintf("Mutual Information: %.4f bits\n", MI_bits))
 
# Example of pre-generated uncorrelated (x1,y1) and correlated variables (x2,y2)
x1 <- as.vector(read.csv("x1.csv",header = FALSE))
y1 <- as.vector(read.csv("y1.csv",header = FALSE))
x2 <- as.vector(read.csv("x2.csv",header = FALSE))
y2 <- as.vector(read.csv("y2.csv",header = FALSE))

source("getMI.R")
MI1 <- getMI(unlist(x1),unlist(y1), nBins=20, smoothingValue = 0.5, units = "bits")
MI2 <- getMI(unlist(x2),unlist(y2),nBins=20, smoothingValue = 0.5, units = "bits")
cat(sprintf("MI of x1 and y1 is  %.3f bits\n", MI1))
cat(sprintf("MI of x2 and y2 is %.3f bits\n", MI2))

