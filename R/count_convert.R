#' A function to extract 'gene_length' from GTF data frame
#'
#' @param gtf_df data.frame formatted GTF file.
#'
#' @importFrom dplyr filter rename
#' @importFrom purrr map
#' @importFrom utils stack
#' @importFrom rlang .data
#'
#' @return feature_lengths in data frame
#' @examples
#' \dontrun{
#'   gtf <- rtracklayer::import('Homo_sapiens.GRCh38.94.gtf')
#'   gtf_df <- as.data.frame(gtf)
#'   gtf_to_gene_length(gtf_df)
#' }
#' @export
gtf_to_gene_length <- function(gtf_df) {
  feature_lengths <- gtf_df %>%
    split(.data$gene_id) %>%
    purrr::map(gtf_to_length) %>%
    stack()
  colnames(feature_lengths) <- c('width', 'ind')
  return(feature_lengths)
}

# helper function
get_gene_length <- function(df) {
  x <- dplyr::filter(df, .data$type == "exon")
  min_start <- min(x$start)
  max_end <- max(x$end)
  bps <- rep(0, max_end - min_start+1)
  for (i in 1:nrow(x))
    bps[x$start[i]:x$end[i] - min_start+1] <- 1
  sum(bps)
}



