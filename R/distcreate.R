#' Create discrete distributions matching given summaries
#'
#' First pass at option 2 at https://github.com/reconhub/covid19hub/issues/4
#'
#' @param name The name of a family of distribution, e.g. `norm`, `gamma`, `weibull`.
#'
#' @param summaries A list of functions calculating distributional summaries
#'   from a sample of the distribution, e.g. `mean`, `sd`, `quantile`.
#'
#' @param values A list of target values to be matched against, each
#'   corresponding to one of the summaries. 
#'
#' @param ini_params A list of named parameters for the distribution, used to
#'   initialise the optimisation procedure.
#'
#' @param weights An optional vector of weights used for the summaries.
#'
#' @param sample_size The size of the sample drawn to compute summaries; larger
#'   values will lead to more accurate estimates, at the cost of memory /
#'   computational time.
#'
#' @param metric An optional metric used to compute dissimilarity between
#'   summaries from the distribution and target values; defaults to the squared
#'   Euclidean distance.
#'
#' @examples
#'
#' if (require(distcrete)) {

#' ini_params <- c(shape = 1, scale = 1)
#' summaries <- c(mean, function(x) quantile(x, c(.25, .75)))
#' values <- list(7, c(4, 11))
#' x <- distcreate("gamma", summaries, values, ini_params = list(shape = 1, scale = 1))
#' summary(x$r(1e5))

#'
#' }


distcreate <- function(name, summaries, values, ini_params, weights = NULL,
                       sample_size = 1e4,
                       metric = function(x, y) (x - y)^2) {

  ## checks
  if (length(summaries) != length(values)) {
    msg <- "`summaries` and `values` have different lengths"
    stop(msg)
  }
  if (!length(summaries)) {
    msg <- "`summaries` is empty"
    stop(msg)
  }


  ## pre-process inputs
  n <- length(summaries)
  values <- unlist(values)
  if (is.null(weights)) {
    weights <- rep(1, length(summaries))
  }
  weights <- unlist(weights)
  
  
  ## make a distcrete object using params
  make_distcrete <- function(params) {
    args <- c(name = name, as.list(ini_params),
              w = 0.5, interval = 1)
    do.call(distcrete::distcrete, args)
  }

  ## get a sample from a distcrete object, return summaries
  get_sample_distance <- function(params) {
    x <- make_distcrete(params)
    samp <- x$r(sample_size)
    out <- rep(0, n)
    for (i in seq_len(n)) {
      out[i] <- sum(metric(summaries[[i]](samp), values[1])) * weights[i]
    }
    sum(out)
  }


  ## run optimisation
  optim_res <- optim(ini_params, get_sample_distance, method = "BFGS")

  ## return results: the distribution, and output of `optim` as attribute
  out <- make_distcrete(optim_res$par)
  attr(out, "optim") <- optim_res
  out
  
}
