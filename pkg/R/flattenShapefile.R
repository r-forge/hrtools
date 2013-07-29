#' Convert 3D Shapefiles into 2D
#'
#' \code{adehabitatHR} works with 2D point locations. Possibly, data in ESRI Shapefile format could contain elevation information, that causes trpuble. This function 'flattens' 3D Shapefiles into a 2D version.
#' A 2D ESRI Shapefile will be written in the same directory (\code{dsn}) containing the input shapefile. 
#' The output shapefile will be named as the input (\code{layer}), plus the "\code{2D}" suffix.
#'
#' @param dsn data source name (the folder containing the 3D ESRI Shapefile)
#' @param layer layer name (the Shapefile name without extension)
#' @return An ESRI Shapefile named as the input (\code{layer}), plus the "\code{2D}" suffix will be written in the same directory (\code{dsn}) of the source spahefile. A \code{SpatialPoint} or \code{SpatialPointDataFrame} is returned.
#' @references \link{http://www.gdal.org/ogr/}, \link{http://www.gdal.org/ogr/ogr_formats.html}
#' @seealso \code{\link{readOGR}} and \code{\link{writeOGR}}
#' @export
#' @examples
#' flat <- flattenShapefile(dsn='~/shapes/', layer='squirrels')
flattenShapefile <- function(dsn=getwd(),layer=NULL) {
  require(rgdal)
  info <- ogrInfo(dsn, layer)
  # is it an ESRI Shapefile (point)?
  stopifnot(info$driver == "ESRI Shapefile", "Input is not an ESRI Shapefile.")
  # is it 3D?
  stopifnot(info$with_z == TRUE, "Input is already 2D")
  # OGRwkbGeometryType is defined as in  ogr_core.h, we need wkbPoint (1), wkbMultiPoint (4), wkbPoint25D (0x80000001), wkbMultiPoint25D (0x80000004)
  stopifnot(info$eType %in% c(1,4), "Input has not a wkbPoint or wkbMultiPoint geometry")
  in.shp <- readOGR(dsn, layer, stringsAsFactors=FALSE, pointDropZ=TRUE)
  out.shpName <- paste(layer,'2D',sep='-')
  writeOGR(in.shp, dsn, out.shpName, driver="ESRI Shapefile", layer_options=c("SHPT=POINT"))
  invisible(in.shp)
}
