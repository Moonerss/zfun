#' @name decompress
#' @title Process multiple type compressed file.
#' @description Process multiple type compressed file. Support zip, tar.gz, tar.bz2, tar, gz
#' @param file Path to compressed file.
#' @param outdir The decompress path.
#' @importFrom R.utils gunzip
#' @importFrom utils untar unzip
#' @export
decompress <- function(file, outdir=NULL) {
  if(is.null(outdir)){
    outdir <- dirname(file)
    cli::cli_alert_info('Decompress file to {outdir}')
  }
  fileType <- basename(file)
  if (grepl("zip$", fileType)) {
    unzip(file, exdir=outdir)
    t1 <- unzip(file, exdir=outdir, list=TRUE)
    return(as.character(t1$Name))
  } else if(grepl("tar\\.gz$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("tar$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("tar\\.bz2$", fileType)) {
    untar(file, exdir=outdir)
    return(untar(file, exdir=outdir, list=TRUE))
  } else if(grepl("\\.gz$", fileType)) {
    gunzip(file, overwrite=TRUE, remove=FALSE)
    if (outdir != dirname(file)) {
      #解压后的文件地址
      t1 <- gsub("[.]gz$", "", file)
      #拷贝到规定的目录
      file.copy(t1, file.path(outdir, basename(t1)))
      file.remove(t1)
      return(basename(t1))
    } else {
      return("OK")
    }
  }
}


#' Read all excel sheets
#'
#' Read specific sheet or all sheets in an excel.
#'
#' @param file Path to the xls/xlsx file.
#' @param sheets name or number to read, if NULL read all sheets.
#' @param merge Whether merge the data, default FALSE.
#' @param verbose Whether print useful message.
#' @param ... Other arguments of \code{\link[readxl]{read_excel}}
#' @importFrom readxl excel_sheets read_excel
#' @importFrom dplyr bind_rows
#'
#' @export
#'
read_all_sheets <- function(file, sheets = NULL, merge = FALSE, verbose = TRUE, ...) {

  if (verbose) cli::cli_alert_info("=> Starting")
  # get sheets
  if (is.null(sheets)) {
    all_sheets <- readxl::excel_sheets(file)
  } else {
    all_sheets <- readxl::excel_sheets(file)
    if (is.numeric(sheets)) {
      all_sheets <- all_sheets[sheets]
    } else if (is.character(sheets)) {
      all_sheets <- intersect(all_sheets, sheets)
    } else {
      cli::cli_abort('`sheets` argument must be a numeric or character vector')
    }
  }

  # read sheets
  if (verbose)
    cli::cli_alert_info("==> Reading sheets: \n", paste0('==> ', paste(all_sheets, collapse = ' ')))
  all_list <- lapply(all_sheets, function(x) {readxl::read_excel(file, sheet = x, ...)})
  names(all_list) <- all_sheets
  if (merge) {
    res <- dplyr::bind_rows(all_list, .id = 'sheet')
  } else {
    res <- all_list
  }
  if (verbose) cli::cli_alert_info("=> Done")
  return(res)
}
