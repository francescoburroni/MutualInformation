#' Calculate mutual information between two variables
#'
#' @description
#' Computes the mutual information I(X;Y) using histogram-based estimation
#' with 10 bins. Smoothing is applied to avoid zero probabilities.
#' Units can be 'nats' (natural log) or 'bits' (log base 2).
#'
#' @param x Numeric vector of observations for variable X
#' @param y Numeric vector of observations for variable Y
#' @param nBins % Positive scalar for bins used in histogram calculations
#' @param smoothingValue Positive scalar for Laplace smoothing to avoid log(0).
#'   Default is 0.5
#' @param units String specifying units: 'nats' or 'bits'. Default is 'nats'
#'
#' @return MI - Mutual information in specified units
#'
#' @examples
#' x <- rnorm(100)
#' y <- x^2
#' MI_nats <- getMI(x, y, smoothingValue = 0.5, units = "nats")
#' MI_bits <- getMI(x, y, smoothingValue = 0.5, units = "bits")
#'
#' @export
getMI <- function(x, y, nBins = 20, smoothingValue = 0.5, units = "bits") {
  
  # Input validation
  # Check that inputs are numeric vectors
  if (!is.numeric(x) || !is.vector(x)) {
    stop("Input 'x' must be a numeric vector")
  }
  if (!is.numeric(y) || !is.vector(y)) {
    stop("Input 'y' must be a numeric vector")
  }
  
  # Check for finite values
  if (any(!is.finite(x))) {
    stop("Input 'x' must contain only finite values")
  }
  if (any(!is.finite(y))) {
    stop("Input 'y' must contain only finite values")
  }
  
  # Check for non-empty vectors
  if (length(x) == 0 || length(y) == 0) {
    stop("Input vectors must not be empty")
  }
  
  # Check smoothingValue
  if (!is.numeric(smoothingValue) || length(smoothingValue) != 1 || 
      smoothingValue <= 0 || !is.finite(smoothingValue)) {
    stop("smoothingValue must be a positive finite scalar")
  }
  
  # Check units
  if (!units %in% c("nats", "bits")) {
    stop("units must be either 'nats' or 'bits'")
  }
  
  # Check that x and y have the same length
  if (length(x) != length(y)) {
    stop("Input vectors x and y must have the same length")
  }
  
  # Check minimum sample size
  if (length(x) < 10) {
    warning(sprintf("Sample size is very small (n=%d). Results may be unreliable.", 
                    length(x)))
  }
  
  # Create bin edges for X and Y spanning their respective ranges
  binEdgesX <- seq(min(x), max(x), length.out = nBins + 1)
  binEdgesY <- seq(min(y), max(y), length.out = nBins + 1)
  
  # Compute marginal histogram for X
  hX <- hist(x, breaks = binEdgesX, plot = FALSE)
  
  # Compute marginal histogram for Y
  hY <- hist(y, breaks = binEdgesY, plot = FALSE)
  
  # Compute joint histogram for (X,Y)
  # Cut data into bins
  x_bins <- cut(x, breaks = binEdgesX, include.lowest = TRUE)
  y_bins <- cut(y, breaks = binEdgesY, include.lowest = TRUE)
  
  # Create frequency table
  hXY_counts <- table(x_bins, y_bins)
  
  # Apply Laplace smoothing to histogram counts to avoid zero probabilities
  XValues <- hX$counts + smoothingValue
  YValues <- hY$counts + smoothingValue
  XYValues <- as.matrix(hXY_counts) + smoothingValue
  
  # Normalize counts to obtain probability distributions
  pX <- XValues / sum(XValues)           # P(X) - marginal distribution
  pY <- YValues / sum(YValues)           # P(Y) - marginal distribution
  pXY <- XYValues / sum(XYValues)        # P(X,Y) - joint distribution
  
  # Create grids for computing independence assumption P(X)P(Y)
  # XGrid[i,j] = P(X=i) for all j
  XGrid <- matrix(rep(pX, each = nBins), nrow = nBins, ncol = nBins, byrow = TRUE)
  
  # YGrid[i,j] = P(Y=j) for all i
  YGrid <- matrix(rep(pY, nBins), nrow = nBins, ncol = nBins, byrow = FALSE)
  
  # Compute product of marginals (independence assumption)
  PInd <- XGrid * YGrid  # PInd[i,j] = P(X=i) * P(Y=j)
  
  # Select logarithm function based on desired units
  if (units == "bits") {
    logFunc <- log2  # Base-2 logarithm for bits
  } else {
    logFunc <- log   # Natural logarithm for nats
  }
  
  # Calculate mutual information using the formula:
  # MI(X;Y) = sum_i sum_j P(X=i,Y=j) * log(P(X=i,Y=j) / (P(X=i)*P(Y=j)))
  MI <- sum(pXY * logFunc(pXY / PInd))
  
  return(MI)
}

# Example usage:
# set.seed(123)
# x <- rnorm(100)
# y <- x^2
# MI_nats <- getMI(x, y, smoothingValue = 0.5, units = "nats")
# MI_bits <- getMI(x, y, smoothingValue = 0.5, units = "bits")
# cat(sprintf("Mutual Information: %.4f nats\n", MI_nats))
# cat(sprintf("Mutual Information: %.4f bits\n", MI_bits))
