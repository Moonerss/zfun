check_args <- function(arg, args) {
  if (!is.null(arg)) {
    arg <- toupper(arg)
    args <- toupper(args)
    if (is.element(arg, args)) {
      return(arg)
    } else {
      cli::cli_abort(c('A bad argument was supplied',
                       'X' = '{.var arg} must be one of {.var {args}}'))
    }
  }
}
