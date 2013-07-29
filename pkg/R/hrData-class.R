#' The Class "HRData": storing multiple Home Ranges
#'
#' A class to store home range calculation results, for multiple home ranges and for multiple estimation methods.
#'
#' The \code{adehabitatHR} package provides several \code{sp}-compatible classes to store home range estimates.
#' 
#' Anyway, when using simultaneously different method (say, \code{mcp} \emph{and} \code{kernelUD}), is is more convenient to have both methods results under the same list element, of course having the reaults organized in a "by animal" list.
#' 
#' The "\code{HRData}" class takes advantage from the "\code{MCHu}" and "\code{estUD}" classes offered by \code{adehabitatHR}, in order to supply a more efficent multiple-method storage class.
#' 
#' @details The \code{HRData} class is basically a list, with an element for each animal, containing all possible results from \code{adehabitatHR} methods.
#' More into detail, a single list element contsins the following slots:
#' \itemize{
#'   \item{\code{geometry}}
#' }
#' 
#' @seealso \link{HRCruncher}, \link[adehabitatHR]{MCHu}, \link[adehabitatHR]{estUD-class}
#
# Some notes on how to structure a data container
# - the idea is that any methos produces, at the end, a SPDF object, from which surface areas can be calculated
# - some particular methods (i.e. kernels) also generate intermediate products such as UDs
# - ?? do we need also some scaffolding for 'statistics"?

# after having skimmed thru https://www.stat.auckland.ac.nz/S-Workshop/Gentleman/S4Objects.pdf and inspird by "Approcah 2" in http://stackoverflow.com/questions/8334796/efficient-way-to-define-a-class-with-multiple-optionally-empty-slots-in-s4-of-r

#### estUDm class is not properly defined in adehabitatHR ######################
setClass("estUDm")

#### HRDataGeom base class, geometries + areas #################################
setClass("HRDataGeom", representation(geometry = "SpatialPolygonsDataFrame", area = "data.frame"))
# union the base class with NULL, to have an instantiable class that can be also NULL
setClassUnion("HRDataGeometry", c("HRDataGeom", "NULL"))

#### HRDataGeom methods
setGeneric("getArea", function(object) standardGeneric("getArea"))
setMethod("getArea", "HRDataGeom", function(object) { 
  a <- data.frame(id=as.character(object@geometry$id), area=object@geometry$area)
  a$id <- as.character(a$id)   # s$id is still factor
  invisible(a)
  }
)

#### HRDataGeomUD, derived class, same as above plus an estUD ##################
setClass("HRDataGeomUD", contains="HRDataGeom", representation(ud = "estUDm"))
setClassUnion("HRDataGeometryUD", c("HRDataGeomUD", "NULL"))

#### HRDataGeomUD methods
setGeneric("geth", function(object) standardGeneric("geth"))
setMethod("geth", "HRDataGeomUD", function(object) {
  # if that UD comes from href, all h are good, if it is a LSCV, return 1 for non-converged h
  res <- sapply(object@ud, function(x) {
    if(x@h$meth == 'LSCV') {
      if(x@h$convergence==TRUE) {
        return(x@h$h)
      } else {
        return(1)
      }
    } else {
      return(x@h$h)
    }
  })
  return(res)
})

#setGeneric("getverticeshr", function(object, ...) standardGeneric("getverticeshr"))
#setMethod("getverticeshr", "HRDataGeomUD", function(object, pct, unin, unout) {
#  object@geometry <- getverticeshr(object@ud, percent=pct, unin=unin, unout=unout)
#})

#### HRData, the final 'superclass' with slots for HRDataGeometry/HRDataGeometryUD objects (or NULL!)
setClass("HRData", representation(mcp="HRDataGeometry", href="HRDataGeometryUD", lscv="HRDataGeometryUD", hadj="HRDataGeometryUD"),
                   prototype(mcp=NULL, href=NULL, lscv=NULL, hadj=NULL))

#@TODO HRDataGeom add method getArea
#@TODO Figure out a method to pull out just the converged lscv UD from an HRDataGeomUD, if it is lscv...
#@TODO HRDataGeomUD add method that does self@geometry <- getverticeshr(self@ud, percent=pct, unin=unin, unout=unout)
