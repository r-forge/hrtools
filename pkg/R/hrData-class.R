#' The Class "HRData": storing multiple Home Ranges
#'
#' A class to store home range calculation results, for multiple home ranges and for multiple estimation methods.
#'
#' The \code{adehabitatHR} package provides several \code{sp}-compatible classes to store home range estimates.
#' 
#' Anyway, when using simultaneously different method (say, \code{mcp} \emph{and} \code{kernelUD}), is is more convenient to all results in a single sata structure, organized as a "by animal ID" list.
#' 
#' The "\code{HRData}" class takes advantage from the "\code{MCHu}" and "\code{estUD}/\code{estUDm}" classes offered by \code{adehabitatHR}, in order to supply a more efficent multiple-method storage class.
#' 
#' @details The \code{HRData} class is an S4 class, containing a slot for esch method. At present, since four methods are supported, the available slots are \code{mcp}, \code{href}, \code{lscv}, \code{hadj}. Each slot is either of "\code{HRDataGeom}" or "\code{HRDataGeomUD}" class. \code{HRDataGeom} is a container for a \code{geometry} slot, which is a \code{SpatialPolygon} or a \code{SpatialPolygonDataFrame}, and an \code{area} slot, a \code{data.frame} containing surface area estimates by animal identifier. \code{HRDataGeomUD} is a derived class with a further \code{ud} slot, used to store utilization distribution data from kernel methods as native \code{estUD}/\code{estUDm} objects.
#' @section Slots:
#' \describe{
#'   \item{\code{mcp}:}{Object of class \code{HRDataGeom}, containing data from \code{mcp} method.}
#'   \item{\code{href}:}{Object of class \code{HRDataGeomUD}, containing data from \code{href} method.}
#'   \item{\code{lscv}:}{Object of class \code{HRDataGeomUD}, containing data from \code{lscv} method.}
#'   \item{\code{hadj}:}{Object of class \code{HRDataGeomUD}, containing data from \code{hadj} method.}
#'  }
#'  @section Methods:
#'  \describe{
#'    \item{\code{getArea(HRData, 'method')}}{returns a \code{data.frame} with home range surface areas.}
#'    \item{\code{getHR(HRData, 'method')}}{returns a \code{SpatialPolygonDataframe} with home range polygons.}
#'    \item{\code{geth(HRData, 'method')}}{returns a \code{numeric} vector with \code{h} values. Note that if \code{method='lscv'}, values of non-converged LSCV estimates will be set to 1.}
#'    \item{\code{hasConverged(HRData)}}{if the \code{lscv} slot is not NULL, returns a logical vector for each animal stating if LSCV mestimator converged or not.}
#'  }
#' @seealso \link{HRCruncher}, \link[adehabitatHR]{MCHu}, \link[adehabitatHR]{estUD-class}
#' @name HRData-class
#
#### estUDm class is not properly defined in adehabitatHR ######################
setClass("estUDm")

#### HRDataGeom base class, geometries + areas #################################
setClass("HRDataGeom", representation(geometry = "SpatialPolygonsDataFrame", area = "data.frame"))
# union the base class with NULL, to have an instantiable class that can be also NULL
## as from https://www.stat.auckland.ac.nz/S-Workshop/Gentleman/S4Objects.pdf,
## and inspired by "Approach 2" in http://stackoverflow.com/questions/8334796/efficient-way-to-define-a-class-with-multiple-optionally-empty-slots-in-s4-of-r
setClassUnion("HRDataGeometry", c("HRDataGeom", "NULL"))

#### HRDataGeom methods
setGeneric("getArea", function(object, ...) standardGeneric("getArea"))
setMethod("getArea", "HRDataGeom", function(object) { 
  a <- data.frame(id=as.character(object@geometry$id), area=object@geometry$area)
  a$id <- as.character(a$id)   # s$id is still factor
  return(a)
  }
)
setGeneric("getHR", function(object, ...) standardGeneric("getHR"))
setMethod("getHR", "HRDataGeom", function(object) {
  return(object@geometry)
})

#### HRDataGeomUD, derived class, same as above plus an estUD ##################
setClass("HRDataGeomUD", contains="HRDataGeom", representation(ud = "estUDm"))
setClassUnion("HRDataGeometryUD", c("HRDataGeomUD", "NULL"))

#### HRDataGeomUD methods
setGeneric("geth", function(object, ...) standardGeneric("geth"))
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
#setMethod("getverticeshr", "HRDataGeomUD", function(object, percent, unin, unout) {
#  object@geometry <- getverticeshr(object@ud, percent, unin, unout)
#})

#### HRData, the final 'superclass' with slots for HRDataGeometry/HRDataGeometryUD objects (or NULL!)
setClass("HRData", representation(mcp="HRDataGeometry", href="HRDataGeometryUD", lscv="HRDataGeometryUD", hadj="HRDataGeometryUD"),
                   prototype(mcp=NULL, href=NULL, lscv=NULL, hadj=NULL))

#' @rdname HRData-class
#' @aliases getArea,HRData,getArea-method
setMethod("getArea", "HRData", function(object, method){
  if(!is.null(slot(object, method))) {
    return(getArea(slot(object, method)))
  }
})

#' @rdname HRData-class
#' @aliases getHR,HRData,getHR-method
setMethod("getHR", "HRData", function(object, method) {
  if(!is.null(slot(object, method))) {
    return(getHR(slot(object, method)))
  }
})

#' @rdname HRData-class
#' @aliases geth,HRData,geth-method
setMethod("geth", "HRData", function(object, method) {
  if(!is.null(slot(object, method))) {
    return(geth(slot(object, method)))
  }
})

#' @rdname HRData-class
#' @aliases hasConverged,HRData,hasConverged-method
setGeneric("hasConverged", function(object, ...) standardGeneric("hasConverged"))
setMethod("hasConverged", "HRData", function(object) {
  if(!is.null(object@lscv)) {
    return(sapply(object@lscv@ud, .hasConverged))
  }
})

