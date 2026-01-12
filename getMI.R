#' Calculate mutual information between two continuous variables
#'
#' @description
#' Computes the mutual information \eqn{MI(X;Y)} using a histogram-based
#' estimator with equally spaced bins and Laplace smoothing applied to the
#' \emph{joint} histogram only. Bin edges are defined as equally spaced
#' intervals between the minimum and maximum of each variable, matching
#' MATLAB's \code{histogram2} behavior. Marginal distributions are derived
#' from the smoothed joint distribution to ensure internal consistency.
#'
#' This estimator is sensitive to the choice of binning and to extreme
#' outliers. In particular, variables with heavy tails may yield low mutual
#' information if most samples fall into a small number of bins.
#'
#' @param x Numeric vector of observations for variable \eqn{X}.
#' @param y Numeric vector of observations for variable \eqn{Y}.
#' @param nBins Positive integer specifying the number of bins per dimension
#'   used in the histogram (default: 20).
#' @param smoothingValue Positive scalar specifying the Laplace smoothing
#'   constant added to each joint histogram cell to avoid zero probabilities
#'   (default: 0.5).
#' @param units Character string specifying the output units. Must be either
#'   \code{"bits"} (base-2 logarithm) or \code{"nats"} (natural logarithm).
#'   Default is \code{"bits"}.
#'
#' @return
#' A single numeric value giving the mutual information between \code{x} and
#' \code{y} in the specified units.
#'
#' @details
#' The mutual information is computed as
#' \deqn{
#' I(X;Y) = \sum_{i,j} p_{ij} \log \frac{p_{ij}}{p_i p_j},
#' }
#' where \eqn{p_{ij}} is the joint probability mass function estimated from
#' a two-dimensional histogram, and \eqn{p_i} and \eqn{p_j} are the marginal
#' distributions obtained by summing the joint distribution.
#'
#' Histogram bins are defined as left-closed, right-open intervals
#' \eqn{[e_i, e_{i+1})}, except for the final bin which includes the right
#' endpoint, reproducing MATLAB's binning convention.
#'
#' @examples
#' set.seed(123)
#'
#' ## Independent variables
#' x1 <- rnorm(1000)
#' y1 <- rnorm(1000)
#' getMI(x1, y1, nBins = 20)
#'
#' ## Nonlinearly dependent variables
#' x2 <- rnorm(1000)
#' y2 <- x2 + 1 / x2
#' getMI(x2, y2, nBins = 20)
#'
#' @export
getMI <- function(x, y, nBins = 20, smoothingValue = 0.5, units = "bits") {

  # ---- Validation ----
  if (!is.numeric(x) || !is.vector(x)) stop("Input 'x' must be a numeric vector")
  if (!is.numeric(y) || !is.vector(y)) stop("Input 'y' must be a numeric vector")
  if (length(x) != length(y)) stop("Input vectors x and y must have the same length")
  if (length(x) == 0) stop("Input vectors must not be empty")
  if (any(!is.finite(x))) stop("Input 'x' must contain only finite values")
  if (any(!is.finite(y))) stop("Input 'y' must contain only finite values")

  if (!is.numeric(nBins) || length(nBins) != 1 || !is.finite(nBins) ||
      nBins != as.integer(nBins) || nBins < 2) {
    stop("nBins must be a single integer >= 2")
  }

  if (!is.numeric(smoothingValue) || length(smoothingValue) != 1 ||
      smoothingValue <= 0 || !is.finite(smoothingValue)) {
    stop("smoothingValue must be a positive finite scalar")
  }

  if (!units %in% c("nats", "bits")) stop("units must be either 'nats' or 'bits'")

  if (length(x) < 10) {
    warning(sprintf("Sample size is very small (n=%d). Results may be unreliable.", length(x)))
  }

  # ---- Constant-vector guard  ----
  if (min(x) == max(x) || min(y) == max(y)) return(0)

  # ---- Bin edges exactly like MATLAB  ----
  edgesX <- seq(min(x), max(x), length.out = nBins + 1)
  edgesY <- seq(min(y), max(y), length.out = nBins + 1)

  # ---- MATLAB-like bin assignment:
  # bins are [edge_i, edge_{i+1}) except last includes right edge
  x_bin <- findInterval(x, edgesX, rightmost.closed = TRUE, all.inside = TRUE)
  y_bin <- findInterval(y, edgesY, rightmost.closed = TRUE, all.inside = TRUE)

  # ---- Joint counts with MATLAB orientation: (Y bins) x (X bins) ----
  XYcounts <- table(
    factor(y_bin, levels = 1:nBins),
    factor(x_bin, levels = 1:nBins)
  )
  XYcounts <- as.matrix(XYcounts)  # nBins x nBins

  # ---- Smooth JOINT only, normalize ----
  pXY <- (XYcounts + smoothingValue)
  pXY <- pXY / sum(pXY)

  # ---- Marginals to match  MATLAB code:
  # pX = sum(pXY,1)' (column sums) ; pY = sum(pXY,2)' (row sums)
  pX <- colSums(pXY)  # length nBins
  pY <- rowSums(pXY)  # length nBins

  # ---- Independence product with same orientation as pXY (Y x X) ----
  PInd <- outer(pY, pX, "*")  # rows=y, cols=x

  logFunc <- if (units == "bits") log2 else log

  MI <- sum(pXY * logFunc(pXY / PInd))

  if (MI < 0 && MI > -1e-12) MI <- 0
  MI
}
