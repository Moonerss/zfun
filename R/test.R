#' Run t test
#' Run t test with two vector after homogeneity of variance test.
#'
#' @param v1 a numeric vector.
#' @param v2 a numeric vector.
#' @param ... Other arguments of \code{\link[stats]{t.test}}
#'
#' @importFrom stats var.test t.test
#'
#' @export
#'
auto_t_test <- function(v1, v2, ...) {
  v1 <- as.numeric(v1)
  v2 <- as.numeric(v2)
  ## Ratio
  if (sum(!is.na(v1)) == 0 | sum(!is.na(v2)) == 0) {
    R_1 <- NA
  } else {
    if (length(v1) == 2 & length(v2) == 2) {
      ## 两次重复的Ratio值
      x <- c(v1, v2)
      a <- c(x[1] / x[3], x[2] / x[4])
      R_1 <- mean(a, na.rm = T)
    } else {
      R_1 <- mean(v1, na.rm = T) / mean(v2, na.rm = T)
    }
  }
  if (sum(!is.na(v1)) < 2 | sum(!is.na(v2)) < 2) {
    P_1 <- NA
  } else {
    cli::cli_alert_info("Do log2 transforming ...")
    v1 <- log2(v1)
    v2 <- log2(v2)
    ## 判断方差齐性
    P_var <- tryCatch(P_var = var.test(v1, v2, alternative = "two.sided")$p.value, error = function(e) {
      P_var <- 1
    })
    if (is.na(P_var)) {
      P_1 <- NA
      # v1 <- c(-9.965784,NA,-9.965784); v2 <- c(-9.965784,-9.965784,NA)  #这种情况P_var=NAN
    } else if (P_var < 0.05) {
      P_1 <- tryCatch(P_1 = t.test(v1, v2, var.equal = F, ...)$p.value, error = function(e) {
        P_1 <- NA
      })
    } else {
      ## p值大于0.05，说明方差齐（没有显著差异）
      P_1 <- tryCatch(P_1 = t.test(v1, v2, var.equal = T, ...)$p.value, error = function(e) {
        P_1 <- NA
      })
    }
  }

  temp <- c(R_1, P_1)

  if (is.na(temp[1])) {
    names(temp) <- c("Ratio", "P_value")
    return(temp)
  } else {
    if (temp[1] < 0.0001) {
      temp[1] <- 0.0001
    }
    temp[1] <- round(temp[1], digits = 3)
    names(temp) <- c("Ratio", "P_value")
    return(temp)
  }
}


#' Run wilcox test
#' Run wilcox test with two vector.
#'
#' @param v1 a numeric vector.
#' @param v2 a numeric vector.
#' @param ... Other arguments of \code{\link[stats]{wilcox.test}}
#'
#' @importFrom stats wilcox.test
#'
#' @export
#'
auto_wilcox_test <- function(v1, v2, ...) {
  v1 <- as.numeric(v1)
  v2 <- as.numeric(v2)

  ## Ratio
  if (sum(!is.na(v1)) == 0 | sum(!is.na(v2)) == 0) {
    R_1 <- NA
  } else {
    if (length(v1) == 2 & length(v2) == 2) {
      ## two repeat ratio
      x <- c(v1, v2)
      a <- c(x[1]/x[3],x[2]/x[4])
      R_1 <- mean(a, na.rm = T)
    } else {
      R_1 <- mean(v1, na.rm = T) / mean(v2, na.rm = T)
    }
  }
  if (sum(!is.na(v1)) < 2 | sum(!is.na(v2)) < 2) {
    P_1 <- NA
  } else {
    message('Do log2 transforming ...')
    v1 <- log2(v1)
    v2 <- log2(v2)
    P_1 <- tryCatch(P_1 = wilcox.test(v1 , v2, ...)$p.value, error = function(e) {P_1 <- NA})
  }
  temp <- c(R_1, P_1)

  if (is.na(temp[1])) {
    names(temp) <- c('Ratio', 'P_value')
    return(temp)
  } else{
    if (temp[1] < 0.0001) {
      temp[1] <- 0.0001
    }
    temp[1] <- round(temp[1], digits = 3)
    names(temp) <- c('Ratio', 'P_value')
    return(temp)
  }
}

# check whether done log2 transformation
# logical value, FALSE not log2, TRUE, log2ed
#' @name check_log2
#' @title Check whether log2 transformed
#' @description \code{check_log2} can check data whether have a log2 transformation
#' @param mat a numeric matrix or a numeric vector for log2 transformation
#'
#' @examples
#' check_log2(1:10)
#' check_log2(matrix(log2(1:10)))
#'
#' @export
#'
check_log2 <- function(mat) {
  qx <- as.numeric(quantile(mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  loged <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (loged) {
    cli::cli_alert_info('The data don\'t log2 transformed!')
    return(FALSE)
  } else {
    cli::cli_alert_info('The data have log2 transformed!')
    return(TRUE)
  }
}
