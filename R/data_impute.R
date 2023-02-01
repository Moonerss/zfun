remove_na <- function(mat, na_ratio = 0.5) {
  cli::cli_alert_info("removing the features with NA more than {na_ratio*100}% precent")
  NA_count <- rowSums(is.na(data))/ncol(data)
  cli::cli_alert_info('remove {sum(NA_count > na_ratio)} features...')
  data <- data[!(NA_count > na_ratio),]
  return(data)
}

fill_by_median <- function(mat) {
  cli::cli_alert_info("Fill the empty value by row median value")
  res <- apply(mat, 1, function(x) {
    x <- as.numeric(x)
    x[is.na(x)] <- median(x, na.rm = T)
    x
  })
  res <- t(res)
  return(res)
}

fill_by_mean <- function(mat) {
  cli::cli_alert_info("Fill the empty value by row mean value")
  res <- apply(mat, 1, function(x) {
    x <- as.numeric(x)
    x[is.na(x)] <- mean(x, na.rm = T)
    x
  })
  res <- t(res)
  return(res)
}

fill_by_min <- function(mat) {
  cli::cli_alert_info("Fill the empty value by the minus value in the data")
  mat[which(is.na(mat))] <- min(mat, na.rm = T)
  return(mat)
}

fill_by_max <- function(mat) {
  cli::cli_alert_info("Fill the empty value by the max value in the data")
  mat[which(is.na(mat))] <- max(mat, na.rm = T)
  return(mat)
}

fill_by_constant <- function(mat, fill_value = 0) {
  cli::cli_alert_info("Fill the empty value by {fill_value}")
  mat[which(is.na(mat))] <- fill_value
  return(mat)
}


#' Fill the NA value by special value
#'
#' Fill the NA value in matrix by special statistical value.
#'
#' @param data A matrix or data.frame contain number.
#' @param fill_type Specific method to used, include:
#' \itemize{
#' \item "median". The NAs will fill by row median.
#' \item "mean". The NAs will fill by row mean.
#' \item "min". The NAs will fill by the minus value in the data.
#' \item "max". The NAs will fill by the max value in the data.
#' }
#' @param fill_value the specific value to fill the NA
#' @param remove Default is FALSE. Whether to remove the genes have too much NA value.
#' @param NA_ratio a number in between 0 and 1, if the NA ratio of a gene greater than it, remove the gene from expression matrix.
#'
#' @importFrom stats median
#'
#' @return A matrix with imputed values
#'
#' @examples
#' \dontrun{
#'   ww <- matrix(data = 1:20, ncol = 4, nrow = 5)
#'   ww[13] <- NA
#'   fill_by_value(data = ww, fill_value = 5)
#'   fill_by_value(data = ww, fill_type = 'min')
#'   fill_by_value(data = ww, fill_type = 'median')
#'   fill_by_value(data = ww, fill_type = 'median', fill_value = 5)
#' }
#' @export
#'
fill_by_value <- function(data, fill_type = NULL, fill_value = NULL,
                          remove = FALSE, NA_ratio = 0.5) {
  ## convert data
  if (is.matrix(data)) {
    data <- data
  } else if (is.data.frame(data)){
    data <- as.matrix(data)
  } else {
    cli::cli_abort(c("x" = '`data` must be matrix or data frame!'))
  }
  ## check argument
  check_args(fill_type, c('mean', 'median', 'min', 'max'))
  ## remove NA genes
  if (remove) {
    data <- remove_na(mat = data, na_ratio = NA_ratio)
  }
  ## fill data
  if (is.null(fill_type)) {
    if (is.null(fill_value)) {
      cli::cli_abort(c('`fill_type` and `fill_value` can\'t all be NULL',
                       "x" = 'You must provide `fill_type` or `fill_value`'))
    } else {
      fill_type <- 'constant'
    }
  } else {
    if (!is.null(fill_value)) {
      cli::cli_alert_warning('`fill_type` and `fill_value` all be provided, use the first one')
    }
  }
  data <- switch(toupper(fill_type),
                 MEAN = fill_by_mean(mat = data),
                 MEDIAN = fill_by_median(mat = data),
                 MAX = fill_by_max(mat = data),
                 MIN = fill_by_min(mat = data),
                 CONSTANT = fill_by_constant(mat = data, fill_value = fill_value)
  )
  return(data)
}

#' Fill NA by normal distribution
#'
#' Fill NA values by random numbers drawn from a normal distribution that has a down-shifted mean
#' and shrunken standard deviation from the sample distribution. This is meant to be similar to imputation
#' in the Perseus software.
#'
#' @param data A matrix or data.frame contain number.
#' @param width Scale factor for the standard deviation of imputed distribution relative to the sample standard deviation.
#' @param downshift Down-shifted the mean of imputed distribution from the sample mean, in units of sample standard deviation.
#' @param remove Default is FALSE. Whether to remove the genes have too much NA value.
#' @param NA_ratio a number in between 0 and 1, if the NA ratio of a gene greater than it, remove the gene from expression matrix.
#' @param log2 Default is FALSE. Whether log2 transformed the data.
#' @param seed Random seed
#'
#' @importFrom stats rnorm sd
#'
#' @return A matrix with imputed values
#'
#' @examples
#' \dontrun{
#'   ww <- matrix(data = 1:20, ncol = 4, nrow = 5)
#'   ww[13] <- NA
#'   fill_by_normal_distribution(ww)
#' }
#'
#' @export
fill_by_normal_distribution <- function(data, width=0.3, downshift=1.8,
                                        remove = FALSE, NA_ratio = 0.5,
                                        log2 = FALSE, seed = 100) {
  ## convert data
  if (is.matrix(data)) {
    data <- data
  } else if (is.data.frame(data)){
    data <- as.matrix(data)
  } else {
    cli::cli_abort(c("x" = '`data` must be matrix or data frame!'))
  }
  ## remove NA genes
  if (remove) {
    data <- remove_na(mat = data, na_ratio = NA_ratio)
  }
  ## check log
  if (log2) {
    data <- log2(data)
  } else {
    mx <- max(data, na.rm=TRUE)
    mn <- min(data, na.rm=TRUE)
    if (mx - mn > 20) {
      cli::cli_alert_warning('Please make sure the values are log-transformed.')
    }
  }
  set.seed(seed)
  data <- apply(data, 1, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  }) %>% t()
  return(data)
}


#' Fill NA by KNN method
#'
#' A function to fill na value using nearest neighbor averaging.
#'
#' @param data A matrix or data.frame contain number.
#' @param remove Default is FALSE. Whether to remove the genes have too much NA value.
#' @param NA_ratio a number in between 0 and 1, if the NA ratio of a gene greater than it, remove the gene from expression matrix.
#' @param ... Other arguments of \code{\link[impute]{impute.knn}}
#'
#' @importFrom impute impute.knn
#'
#' @return A matrix with imputed values
#'
#' @examples
#' \dontrun{
#'   ww <- matrix(data = 1:20, ncol = 4, nrow = 5)
#'   ww[13] <- NA
#'   fill_by_knn(ww)
#' }
#'
#' @export
fill_by_knn <- function(data, remove = FALSE, NA_ratio = 0.5, ...) {
  ## convert data
  if (is.matrix(data)) {
    data <- data
  } else if (is.data.frame(data)){
    data <- as.matrix(data)
  } else {
    cli::cli_abort(c("x" = '`data` must be matrix or data frame!'))
  }
  ## remove NA genes
  if (remove) {
    data <- remove_na(mat = data, na_ratio = NA_ratio)
  }
  res <- impute.knn(data, ...)$data
  colnames(res) <- colnames(data)
  return(res)
}
