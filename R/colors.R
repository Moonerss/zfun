#' @title Collected general color palette
#'
#' @description Many beautiful color palette collected by Erjie.
#'
#' @param n number of colors to generate
#' @param alpha factor modifying the opacity alpha; typically in [0,1]
#' @param theme color from different color theme
#'
#' @importFrom grDevices adjustcolor
#' @importFrom RColorBrewer brewer.pal
#'
#' @return return a color vector
#' @examples
#' get_color_palette(n = 5, alpha = 0.7, theme = 'cellchat')
#'
#' @export
#'
get_color_palette <- function(n, alpha = 1, theme = c('maftools_domain_col', 'maftools_mutation1',
                                                      'maftools_mutation2', 'maftools_snv',
                                                      'cellchat', 'palette1')) {
  ## get colors from different articles or packages
  theme <- match.arg(theme)
  ## from maftools package
  maftools_mutation_col1 <- c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
                              "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
                              "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
                              "#607D8B")
  maftools_mutation_col2 <- c(RColorBrewer::brewer.pal(11, name = "Paired"),
                              RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],
                              'black', 'violet', 'royalblue', '#7b7060', '#535c68')
  maftools_snv_col <- c("#F44336", "#3F51B5", "#2196F3", "#4CAF50", "#FFC107", "#FF9800")
  maftools_domain_col <- c("#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#f19066",
                           "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f8a5c2",
                           "#63cdda", "#ea8685", "#596275", "#574b90", "#f78fb3", "#3dc1d3",
                           "#e66767", "#303952")
  # palette1 <- c(brewer.pal(8,"Set1")[-6], rev(brewer.pal(8,"Dark2")), brewer.pal(7,"Set2"))[c(1:12,16:19,13:15)]
  palette1_col <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
                "#F781BF", "#666666", "#A6761D", "#E6AB02", "#66A61E", "#E7298A",
                "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#7570B3", "#D95F02", "#1B9E77")
  ## from cellchat
  cellchat_col <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF',
                    '#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                    '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3',
                    '#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69',
                    '#8DD3C7','#999999')

  ## from article
  # 1. a simple article

  col <- switch (theme,
                 maftools_domain_col = maftools_domain_col,
                 maftools_mutation1 = maftools_mutation_col1,
                 maftools_mutation2 = maftools_mutation_col2,
                 maftools_snv = maftools_snv_col,
                 cellchat = cellchat_col,
                 palette1 = palette1_col
  )
  ## set alpha
  col <- grDevices::adjustcolor(col = col, alpha.f = alpha)

  ## return
  if (n <= length(col)) {
    colors <- col[1:n]
  } else {
    warning('The ', theme, ' theme only contain ', length(col), ' colors, we keep all.')
    colors <- col
  }
  return(colors)
}

