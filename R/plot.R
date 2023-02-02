#' @title Make Venn Diagram from a list
#' @description Make Venn Diagram base on \link[VennDiagram:venn.diagram]{venn.diagram}, you can save it by \code{ggsave}
#' @param x A list of vectors (e.g., integers, chars), with each component corresponding to a separate circle in the Venn diagram
#' @param ... other argument of \link[VennDiagram:venn.diagram]{venn.diagram}
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid textGrob gList grid.newpage grid.draw
#' @export
#' @author Erjie Zhao
#' @seealso \link[VennDiagram:venn.diagram]{venn.diagram}
#' @examples
#' \dontrun{
#'   kp <- plot_venn(list(a = sample(letters[1:10], 20, replace = T),
#'                        b = sample(letters[1:10], 25, replace = T)),
#'                   fill=c("#CC79A7", "#56B4E9"),
#'                   col=c("#D55E00", "#0072B2"),
#'                   cat.col=c("#D55E00", "#0072B2"))
#' }
#' @export
plot_venn <- function(x, ...) {
  ## set base attributes
  dots <- list(...)
  if(is.null(dots$cat.cex)) dots$cat.cex = 1
  if(is.null(dots$cat.col)) dots$cat.col = "black"
  if(is.null(dots$cat.fontface)) dots$cat.fontface = "plain"
  if(is.null(dots$cat.fontfamily)) dots$cat.fontfamily = "serif"

  dots$x <- x
  dots <- c(dots, filename = list(NULL))
  venngrid <- do.call(VennDiagram::venn.diagram, dots)
  unlink(dir(pattern="^VennDiagram.[0-9_\\-]+.log$"))
  ## show picture in windows
  if (Sys.info()['sysname'] == "Windows") {
    grid.newpage()
    grid.draw(venngrid)
  }
  return(venngrid)
}
