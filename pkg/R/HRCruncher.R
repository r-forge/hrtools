#' Batch home range calculations.
#'
#' Driver function for batch homerenge calculations.
#'
#' @param xy An object inheriting the class \code{SpatialPoints} containing the x and y relocations of the animal. If \code{xy} inherits the class \code{SpatialPointsDataFrame}, it should contain a column whose values specify the identity of the animals for each relocation. If that column exists and is not named \code{ID}, its name must be specified using the \code{id} argument. If a \code{SpatialPoints} object is passed, it is assumed that all the relocations are from a single animal.
#' @param method a string specifying the home range calculation method to be used. See Details for further specifications.
#' @param grid a number giving the \emph{cell size} (in meters) of the grid on which the UD should be estimated. Default is 100 m. Alternatively, this parameter may be an object inheriting the class SpatialPixels, that can be created using the \code{sizeGrid} function.
#' @param idfield the name of a column in \code{xy} corresponding to the identity of the animals for each relocation.
#' @param minfix an integer specifying the minimum number of relocations needed to calculate an home range. Default value is 15. If an animal has less than \code{minfix} relocations a warning will be issued and no home renge calculations will be made for that animal.
#' @param percent for \code{mcp} method, the number of relocations to use to calculate a Minimum Convex Polygon. The default value of 95 means that a 5% of points falling far from the centroid of the relocations used will be treated as outliers. See \code{percent} argument of \code{mcp} function. For other (kernel-based) methods \code{prc} represents the percentage level (default: 95%) for home range estimation (see \code{percent} parameter in \code{kernelUD} function).
#' @param unin the units of the relocations coordinates. Either "m" for meters (default) or "km" for kilometers.
#' @param unout the units of the output areas. Either "m2" for square meters, "km2" for square kilometers or "ha" for hectares (default).
#' @param boundary If not NULL, an object inheriting the class SpatialLines defining a barrier that cannot be crossed by the animals. See \code{kernelUD} for details.
#' @return an object of \code{HRData} class, with its slots filled.
#' @details Available methods (see \code{method} above) are: \code{mcp}, \code{href}, \code{lscv}, \code{hadj}.
#' \itemize{
#'   \item{\code{mcp}}{ Minimum Convex Polygon, using \code{adehabitatHR} \code{mcp} function.}
#'   \item{\code{href}}{ kernel home range using \code{adehabitatHR} \code{kernelUD} function with "ad hoc" smoothing parameter \code{h}.}
#'   \item{\code{lscv}}{ kernel home range using \code{adehabitatHR} \code{kernelUD} function with least-squares cross validation to estimate \code{h} smoothing parameter.}
#'   \item{\code{hadj}}{ kernel home range using \code{adehabitatHR} \code{kernelUD} function with "adjusted" \code{h} as in Wauters et al. (2007).}
#'   }
#' @seealso \link{sizeGrid}, \link[adehabitatHR]{mcp}, \link[adehabitatHR]{kernelUD}, \link{HRData-class}
#' @references Lucas A. Wauters, Damiano G. Preatoni, Ambrogio Molinari, Guido Tosi (2007).
#' Radio-tracking squirrels: Performance of home range density and linkage estimators with small range and sample size.
#' Ecological Modelling 202(3-4): 333-344 \url{http://dx.doi.org/10.1016/j.ecolmodel.2006.11.001}
#' @examples
#' data(squirrels)
#' # size up a grid with 30 m cell size
#' reference.grid <- sizeGrid(squirrels, res=30)
#' 
#' # plot reference grid and locations
#' image(reference.grid, col='lightgray')
#' plot(squirrels, col=seq(1:18), add=TRUE)
#' 
#' # define animal ID field
#' idfield <- "TAG"
#' 
#' # define minumum number of fixes needed
#' minfix <- 15
#' 
#' calculate mcp and hadj kernel (href and lscv kernels will be calculated implicitly)
#' result <- HRCruncher(squirrels, method=c('mcp','hadj'), percent=95, grid=reference.grid, idfield=idfield, unin='m', unout='ha')
#' 
#' # save results as shapefile
#' libaray(rgdal)
#' writeOGR(result@@href@@geometry, dsn='/tmp', layer='test' , driver="ESRI Shapefile")
HRCruncher <- function(xy, method=c("mcp","href","lscv","hadj"), grid=100, idfield=c("ID"), minfix=15, percent=95, unin=c("m","km"), unout=c("ha", "km2", "m2"), boundary=NULL) {
  # methods available in adehabitatHR, by input data type, in alphabetical order
  # input is a SpatialPoints or a SpatialPointsDataFrame:
  # CharHull(xy - sp)
  # clusthr(xy - sp)
  # kernelUD(xy - sp) i.e. href, lscv, hadj
  # LoCoH.k(xy - sp)
  # LoCoH.r(xy - sp)
  # LoCoH.a(xy - sp)
  # mcp(xy - sp)
  # input is a ltraj (trajectory)
  # BRB(ltraj)
  # BRB.D(ltraj)
  # BRB.likD(ltraj)
  # kernelbb(ltraj)
  # kernelkc(ltraj)
  # input is a data frame with x,y,time
  # kernelkcbase(xyt)
  # here we limit to 'classic' xy-based methodologies, i.e. non-trajectory based methods.
  # mcp, href, lscv, hadj, clusthr, CharHull, LoCoH.k, LoCoH.r, LoCoH.a
  
  # check here for method, grid, id, minfix, prc, unin, unout using either missing() or match.arg()
  # actually, we chack only for xy presence and for valid values for method, unin, unout
  if(missing(xy)) {error("xy argument is mandatory.")}
  method <- match.arg(method, several.ok = TRUE)
  unin <- match.arg(unin)
  unout <- match.arg(unout)
  
  # check what xy actually is
  if(inherits(xy, "SpatialPoints")) {
    # xy either is a SpatialPoints or a SpatialPointsDataFrame
    if(!inherits(xy, "SpatialPointsDataFrame")) {
      # xy _is_ a SpatialPoints, just add the DataFrame part with a constant ID field
      data <- data.frame(ID = character(length(xy)))
      data$ID <- '0'
      xy <- SpatialPointsDataFrame(xy, data)
    }
    # now we're sure that xy is a SpatialPointsDataFrame, with _at least_ an ID field.
    # actually we need a SpatialPointsDataFrame with _exactly one_ field, called ID: check that and fix if needed
    if(ncol(xy) == 1) {
      # data slot has exactly 1 column: assume that is ID
      names(xy) <- 'ID'
    } else {
      # data slot has more than one column: look for an user-specified ID field
      if (idfield %in% names(xy)){
        # store attributes in a separate dataframe useful to append to hr at the end of the process
        #@TODO mind that 'attributes' can be re-attached via ID field only. this means that is useless to save 'attributes' into dd as they are in the SPDF, we must pass them thru unique(), unicizing them by ID. Inecitably, if there are different value sets per ID repetitions, we will lose them...
        dd <- xy@data
        # rearrange to 'one column SPDF'
        xy <- subset(xy, select=(idfield))
        names(xy) <- 'ID'
      } else {
        stop("a field named '", idfield, "' was not found. Please check input data.")
      }
    }
  } else {
    #@TODO check whether we have an 'old style' tabular structure in hand, and convert it, place code here. At present, we just issue an error and abort.
    #@NOTE(Fra) adehabitatHR gives 2 fucntions for Conversion of old classes from adehabitat to classes from adehabitatHR (kver2spol; khr2estUDm)
    #if (class(x) == "sahrlocs") {
    #  dd <- getsahrlocs(x, what="locs")
    #  dd$ID <- dd$Name
    #  dd$Name <- NULL
    #}
    #@NOTE(prea) these seem to be functions to convert a _result_. we instead need to convert a dataframe with x, y, ID into a SPDF.
    stop("xy should be of class SpatialPoints or SpatialPointsDataFrame")
  }
  # at this point we have a SpatialPointsDataFrame with a single 'ID' field in its data slot.
  # to be sure, force ID to be of type character (character, _not_ factor)
  xy$ID <- as.character(xy$ID)
  
  # subset animals that have at least minfix fixes
  counts <- as.data.frame(xtabs(~ ID, xy))
  ID_to_remove <- as.character(subset(counts, Freq < minfix)[,"ID"])
  if(length(ID_to_remove > 0)) {
    xy <- xy[!xy$ID %in% ID_to_remove,]
    warning("The following animals have been excluded from home range calculations, having less than ", minfix, " relocations:\n\t", paste(ID_to_remove, collapse=', '), ".")
  } 
  
  # actual hr calculations here ################################################
  # note: only geometries are calculated, that is, for each method we end up 
  # with a SpatialPolygonsDataFrame, and any intermadiate product (e.g. UDs).
  # so, any method will for sure yield a <method>GeomRes (SPDF) and some methods
  # will also yield a <method>UDRes.
  # there is no need to calculate SPDF surface areas, since SPDF already carries
  # them in data and id slots. an internal function .getArea(SPDF) has been 
  # devised to return a dataframe with id and area
  
  # create a container for results. See hrData-class.R, HRRes will have all its slots set to NULL
  HRRes <- new("HRData")

  # mcp
  if ("mcp" %in% method) {                                                ## MCP
    cat("MCP calculation running...\n")
    #mcpGeomRes <- mcp(xy, percent=pct, unin=unin, unout=unout)
    HRRes@mcp <- new("HRDataGeom", geometry=mcp(xy, percent=percent, unin=unin, unout=unout))
    # to calculate surface area, calling mcp.area is useless, singe we already gor geometries. 
    # Not Run: mcpAreaRes <- .getArea(mcpGeomRes)
    # BETTER: the HRDataGeom base class will have a getArea method...
    cat("Done.\n")  
  }
  
  # href
  if ("href" %in% method | "hadj" %in% method) {                    ## KHR(href)
    cat("Kernel href calculation running...\n")
    HRRes@href <- new("HRDataGeomUD", ud=kernelUD(xy, h="href", grid=grid, same4all=FALSE, kern="bivnorm", boundary=boundary))
    HRRes@href@geometry <- getverticeshr(HRRes@href@ud, percent=percent, unin=unin, unout=unout)
    cat("Done.\n")
  } 
  
  # lscv
  if ("lscv" %in% method | "hadj" %in% method) {                   ## KHR(lscv)
    cat("Kernel LSCV calculation running...\n")
    options(warn=-1) # turn off warnings, since LSCV can output lots of almost useless stuff when not converging
    #lscvUDRes <- kernelUD(xy, h="LSCV", grid=grid, same4all=FALSE, kern="bivnorm", boundary=boundary) 
    HRRes@lscv <- new("HRDataGeomUD", ud=kernelUD(xy, h="LSCV", grid=grid, same4all=FALSE, kern="bivnorm", boundary=boundary))
    options(warn=0) # be warned again...
    # for some UD, LSCV could have not converged, and getverticeshr would return errors:? check lscvUDRes@h$convergence
    nc <- sapply(HRRes@lscv@ud, .hasConverged)
    if(!all(nc)) {
      # at least one LSCV did not converge: filter out and warn in a more informative way than kernelUD
      warning("The least-square cross-validation algorithm did not converge for the following animals:\n\t", paste(names(nc[nc==FALSE]), collapse=', '),"\nleast-square cross-validated kernels will not be calculated for those animals. ")
      lscvUDRes <- Filter(.hasConverged, HRRes@lscv@ud)
      class(lscvUDRes) <- "estUDm" # Filter changes lscvUDRes class from estUDm to List
    }
    HRRes@lscv@geometry <- getverticeshr(lscvUDRes, percent=percent, unin=unin, unout=unout)
    cat("Done.\n")
  }
  
  # hadj 
  if ("hadj" %in% method) {                                          ## KHR(hadj)
    cat("Kernel hadj calculation running...\n")
    # hadj calculation is a three-step process: 
    #  1) calculate (where possible) hadj[i] = hlscv[i]/href[i], for each animal i
    #  2) calculate average hadj.mean 
    #  3) recalculate fixed h kernel surface areas with h[i] = hadj.mean * href[i]
    hadji <- geth(HRRes@lscv) / geth(HRRes@href)
    hadjm <- mean(hadji)
    hadj <- hadjm * geth(HRRes@href)
    # kernelUD (as well as other adehabitatHR internals) does not accept a vector of h values
    # a solution has been proposed in http://lists.faunalia.it/pipermail/animov/2013-April/001073.html
    xyList <- split.data.frame(xy, xy$ID) # call the right 'split' method 
    hadjUDRes <- lapply(seq_along(xyList), function(x) kernelUD(xyList[[x]], hadj[x], grid=grid, same4all=FALSE, kern="bivnorm", boundary=boundary))
    hadjUDRes <- lapply(hadjUDRes, function(x) x[[1]]) # dirty fix because we got a (by-id) list of lists with just one element...
    names(hadjUDRes) <- names(xyList) # because the dirty fix destroyed names
    class(hadjUDRes) <- "estUDm"
    HRRes@hadj <- new("HRDataGeomUD", ud=hadjUDRes)
    HRRes@hadj@geometry <- getverticeshr(HRRes@hadj@ud, percent=percent, unin=unin, unout=unout)
    # Not Run: hrefAreaRes <- .getArea(hrefGeomRes)
    cat("Done.\n")
  }
    
  invisible(HRRes)
} ## end of HRCruncher function

## used to check convergence achievement in LSCV kernel UDs, pplies to estUD object (NOT estUDm!)
.hasConverged <- function(x) {
  return(x@h$convergence==TRUE)
}
