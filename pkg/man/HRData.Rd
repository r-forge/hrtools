\name{HRData}
\alias{HRData}
\title{The Class "HRData": storing multiple Home Ranges}
\description{
  A class to store home range calculation results, for
  multiple home ranges and for multiple estimation methods.
}
\details{
  The \code{adehabitatHR} package provides several
  \code{sp}-compatible classes to store home range
  estimates.

  Anyway, when using simultaneously different method (say,
  \code{mcp} \emph{and} \code{kernelUD}), is is more
  convenient to have both methods results under the same
  list element, of course having the reaults organized in a
  "by animal" list.

  The "\code{HRData}" class takes advantage from the
  "\code{MCHu}" and "\code{estUD}" classes offered by
  \code{adehabitatHR}, in order to supply a more efficent
  multiple-method storage class.

  The \code{HRData} class is basically a list, with an
  element for each animal, containing all possible results
  from \code{adehabitatHR} methods. More into detail, a
  single list element contsins the following slots:
  \itemize{ \item{\code{geometry}} }
}
\seealso{
  \link{HRCruncher}, \link[adehabitatHR]{MCHu},
  \link[adehabitatHR]{estUD-class}
}

