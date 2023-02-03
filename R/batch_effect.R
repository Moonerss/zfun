#' @title Remove batch effect of RNA-seq data
#'
#' @description Remove batch effect of RNA-seq data by combat
#'
#' @param expressionMatrix Genomic measure matrix (dimensions probe x sample) - for example, expression matrix.
#'        It is better to use a log transformed data as input, and there must be no negative value in the data,
#'        because it will generate NaN value.
#' @param batch Batch covariate (only one batch allowed)
#' @param covariate A data frame contain the covariates instead in. Each column must be factor with the same sample
#'        order with the data. For example, If the batch of different experimental platforms is removed and the
#'        differences of normal and cancer are calculated, then the normal and a CANCER samples in each experiment
#'        are variables that are interested, this data are generally variables that are interested in subsequent
#'        research. When batch correction, you need to consider the difference in quantity distribution in each batch (unbalanced)
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param ... Other arguments of \code{\link[sva]{ComBat}}
#' @importFrom stats model.matrix
#'
#' @seealso
#' \url{https://www.biostars.org/p/196430/}
#' \url{https://support.bioconductor.org/p/73266/}
#' \url{https://www.biostars.org/p/173684/}
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#'
#' @note If the result of \code{Combat} contain negative values, replace by zero
#'
#' @examples
#' \dontrun{
#'
#' library(bladderbatch)
#' data(bladderdata)
#' pheno = pData(bladderEset)
#' expressionMatrix = exprs(bladderEset)
#' batch = pheno$batch
#' covariate = data.frame(cancer = as.factor(pheno$cancer))
#' correctMatrix <- FunctionCombat(expressionMatrix, batch, covariate = covariate)
#'
#' }
#' @export
#'
batch_combat <- function(expressionMatrix, batch, covariate = NULL, par.prior = TRUE, prior.plots = FALSE, ...){

  # 由于当前只进行batch的矫正，所以建立null model
  phen <- data.frame(sample = 1:ncol(expressionMatrix), batch = batch)
  rownames(phen) <- colnames(expressionMatrix)

  if(is.null(covariate)){
    mod <- model.matrix(~1, data = phen)
  }else{
    dataType <- sapply(covariate, as.factor)
    mod <- model.matrix(~., data = covariate)
  }

  combat_edata <- sva::ComBat(dat=expressionMatrix, batch=batch, mod=mod, par.prior=par.prior, prior.plots = prior.plots, ...)

  # 将矫正后的表达值为负数的设置为0
  if (any(combat_edata<0)) {
    cli::cli_alert_warning('There have nagative value, replace by 0.')
  }
  combat_edata[which(combat_edata<0)] <- 0

  return(combat_edata)
}
