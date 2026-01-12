# Example usage:
set.seed(1310)
x <- rnorm(1000)
y <- x
source("getMI.R")
MI_nats <- getMI(x, y, nBins=20, smoothingValue = 0.5, units = "nats")
MI_bits <- getMI(x, y, nBins=20, smoothingValue = 0.5, units = "bits")
cat(sprintf("Mutual Information: %.4f nats\n", MI_nats))
cat(sprintf("Mutual Information: %.4f bits\n", MI_bits))
 