# median_mad
median_mad <- function(matrix) {
  norm.params <- t(apply(matrix, 2, function (x) c(center = median(x, na.rm = T), scale = mad(x, na.rm = T))))
  data.norm <- scale(matrix, center = norm.params[,'center'], scale = norm.params[,'scale'])
  return(data.norm)
}

# mean_sd
mean_sd <- function(matrix) {
  norm.params <- t(apply(matrix, 2, function (x) c(center = mean(x, na.rm = T), scale = sd(x, na.rm = T))))
  data.norm <- scale(matrix, center = norm.params[,'center'], scale = norm.params[,'scale'])
  return(data.norm)
}

# quantile
quantile_normalize <- function(matrix) {
  data.norm <- normalize.quantiles(matrix)
  return(data.norm)
}

# upper quantile
upper_quantile_normalize <- function(matrix) {
  data.norm <- apply(matrix, 2, function(x) x - quantile(x, c(0.75),na.rm=T))
  return(data.norm)
}

# VSN
vsn_normalize <- function(matrix) {
  data.norm <- justvsn(matrix)
  return(data.norm)
}

# two-component mixture model normalization
two_comp <- function (sample, type = c('default', 'bimodal'), mode.lower.bound=-3) {
  #   1. For all sample types, fit a 2-component gaussian mixture model using normalmixEM.
  #   2. For the bimodal samples, find the major mode M1 by kernel density estimation
  #     2a. Fit the model with one component mean constrained to equal M1
  #     2b. Normalize (standardize) samples using mean (M1) and resulting std. dev.
  #   3. For unimodal (default) samples, find the mode M using kernel density estimation
  #     3a. Fit the model with mean for both components constrained to be equal to M
  #     3b. Normalize (standardize) samples using mean M and smaller std. dev. from model fit
  #
  #  the major mode should be located at a value larger than mode.lower.bound

  # WARNING:
  # This code has a lot of hacks to fix the flakiness of normalmixEM, and the idiosyncracies
  # of the actual data. Carefully re-examine code for new or altered input data
  # Currently works for log-ratio data approximately centered around 0
  is.error <- function(x) inherits(x, "try-error")             # function to check for error
  type <- match.arg(type)

  data <- sample [ !is.na (sample) ]
  data.range <- diff (range (data))
  dens <- try (density (data, kernel='gaussian', bw='SJ'))     # gaussian kernel with S-J bandwidth
  if (is.error (dens))                                         # sometimes, SJ bw estimation fails
    dens <- density (data, kernel='gaussian', bw='ucv')        # in such cases, use unbiased CV
  # (see Venalbles & Ripley, 2002, pg, 129
  #  and Density Estimation, S.J.Sheather, Stat. Sci. 2004)
  # find major (highest) mode > -3 (to avoid problems with lower mode having higher density than higher mode)
  x.range <- dens$x > mode.lower.bound
  dens.x <- dens$x [x.range];  dens.y <- dens$y [x.range]
  mode <- dens.x[which.max(dens.y)]
  if (type=='bimodal') mean.constr <- c (NA, mode) else mean.constr <- c (mode, mode)
  model <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
  model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
  model.alt <- Mclust (data, G=2, modelNames=c ("V","E"))
  # V results is separate SDs for each cluster; E fits a single SD for both clusters
  if (length (model.alt$parameters$variance$sigmasq)==1)  # latter code expects two SD values
    model.alt$parameters$variance$sigmasq <- rep (model.alt$parameters$variance$sigmasq, 2)
  alt.mu <- model.alt$parameters$mean
  alt.sd <- sqrt (model.alt$parameters$variance$sigmasq)
  # find reproducible model fit that is close to Mclust fit
  # if not, re-fit model -- without this condition
  # normalmixEM produces one-off model fits
  n.try <- 1
  if (type=='bimodal') model.mode <- which(model$mu==mode)
  else model.mode <- which(model$mu==mode)[which.min (model$sigma)]
  model.other <- model.mode %% 2 + 1
  alt.mode <- ifelse (diff (alt.mu) < data.range*0.05,          # means are close --
                      which.min (alt.sd),                       # use min sd to pick alt.mode
                      which.min(abs(model.alt$par$mean-mode)))  #  else use alt.mu closest to mode
  # always using latter can result in consistently picking the wrong alt.mode => no convergence
  alt.other <- alt.mode %% 2 + 1
  while ( abs (model$mu[model.mode] - alt.mu[alt.mode]) > data.range*0.05 ||
          abs (model$sigma[model.mode]-alt.sd[alt.mode]) > data.range*0.05 ||
          model$sigma[model.mode] < min (0.1, data.range*0.01) ||   # some samples can have very small SD (eg CR replicates)
          (type=='bimodal' && (abs (model$mu[model.other] - alt.mu[alt.other]) > data.range*0.25)) ||
          abs (sum (c (model$mu, model$sigma) - c (model.rep$mu, model.rep$sigma))) > 1e-3 ) {
    # if major mode (and SD of mode) is not within 5% of data range, or if the other mean (for bimodals only)
    # is not within 25% of the Mclust result, try again
    model <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)
    model.rep <- normalmixEM (data, k=2, mean.constr=mean.constr, maxit=10000)

    if (n.try > 50) stop (paste ("Can't fit mixture model ... giving up\n"))
    n.try <- n.try + 1
  }


  if (type=='bimodal') {
    # sometimes (esp. in phosphoproteome) the minor (lower) mode can be larger than the major (higher) mode
    # this situation is not possible in the unimodal samples
    corrected.mode <- model$mu [which.max(model$mu)]
    if (corrected.mode != mode) {
      cat ('  Lower mode larger than higher mode\n')
      mode <- corrected.mode
    }
  }
  norm.mean <- mode
  norm.sd <- ifelse (type=='bimodal', model$sigma[which(model$mu==mode)], min (model$sigma))

  # normalize by standardizing
  data <- data - norm.mean
  data <- data / norm.sd

  # return normalized data reorganized to original order
  sample [ !is.na (sample) ] <- data
  return ( list (norm.sample=sample, norm.mean=norm.mean, norm.sd=norm.sd, fit=unlist (c(model$mu, model$sigma))) )
}

two_comp_normalize <- function(matrix) {
  data.norm.list <- vector('list', ncol(matrix))
  names(data.norm.list) <- colnames(matrix)
  for(x in colnames(matrix)){
    res <- try(two_comp(matrix[, x], type="default"))
    data.norm.list[[x]] <- res
    if(inherits(res, 'try-error')) break;
  }
  ## if 2-comp was successful on all data column convert list to matrix
  data.norm <- matrix(unlist(lapply(data.norm.list, function(x)x$norm.sample)),
                      ncol=length(data.norm.list),
                      dimnames=list(rownames(matrix), names(data.norm.list)))
  return(data.norm)
}

#' Normalize the expression data
#'
#' Normalize expression data by different scale method.
#'
#' @param matrix The expression matrix with gene in row and sample in column.
#' @param method The method used to scale data, include:
#' \itemize{
#' \item "median_mad". Normalize data by .
#' \item "mean_sd". The NAs will fill by row mean.
#' \item "quantile". Normalize data by quantile method in \code{\link[preprocessCore]{normalize.quantiles}}.
#' \item "upper_quantile". Normalize data by upper 75 percent quantile.
#' \item "vsn". Normalize data by variance stabilizing normalization method in \code{\link[vsn]{justvsn}}.
#' \item "2-component". Normalize data by two-component mixture model of Michael A. Gillette et.al.
#' }
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom vsn justvsn
#' @importFrom mixtools normalmixEM
#' @importFrom mclust Mclust
#' @importFrom stats density mad quantile
#'
#' @references
#' [1] Proteogenomic Characterization Reveals Therapeutic Vulnerabilities in Lung Adenocarcinoma.
#'
#' @return return a normalized matrix
#' @examples
#' \dontrun{
#'   data(kidney_expr)
#'   normalize_data(kidney, method = 'vsn')
#' }
#' @export
#'

normalize_data <- function(matrix, method = c('median_mad', 'mean_sd', 'quantile', 'upper_quantile', 'vsn', '2-component')) {
  method <- match.arg(method)
  cli::cli_alert_info('Normalize method: {method}')
  data_norm <- switch(method,
                      median_mad = median_mad(matrix),
                      mean_sd = mean_sd(matrix),
                      quantile = quantile_normalize(matrix),
                      upper_quantile = upper_quantile_normalize(matrix),
                      vsn = vsn_normalize(matrix),
                      `2-component` = two_comp_normalize(matrix)
                      )
  rownames(data_norm) <- rownames(matrix)
  colnames(data_norm) <- colnames(matrix)
  return(data_norm)
}

