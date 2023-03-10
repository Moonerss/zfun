% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verse.R
\name{zfun_conflicts}
\alias{zfun_conflicts}
\title{Conflicts between the zfun and other packages}
\usage{
zfun_conflicts()
}
\description{
This function lists all the conflicts between packages in the zfun
and other packages that you have loaded.
}
\details{
If dplyr is one of the select packages, then the following four conflicts
are deliberately ignored: \code{intersect}, \code{union}, \code{setequal},
and \code{setdiff} from dplyr. These functions make the base equivalents
generic, so shouldn't negatively affect any existing code.
}
\examples{
zfun_conflicts()
}
