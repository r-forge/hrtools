#' Create a reference grid for kernel calculations
#'
#' Create a \code{SpatialPixels} object with same extent and CRS of input points and resolution as desired.
#'
#' @param xy an object inheriting the class \code{SpatialPoints} containing the x and y relocations of the animal. That object will be used to calculate grid extent and CRS.
#' @param res The grid spatial resolution, in meters.
#' @param as specify whether to output the reference grid as \code{SpatialPixels} object (default) or as \code{raster}.
#' @return a \code{SpatialPixels} object (default) or a \code{raster} object.
#' @examples
#' # generate a reference grid with a 50 m spatial resolution for the squirrel dataset
#' data(squirrels)
#' reference.grid <- sizeGrid(squirrels, res=50)
sizeGrid <- function(xy, res=100, as=c("SpatialPixels","raster")) {
  if (!inherits(xy, "Spatial")) 
    stop("xy should inherit the class Spatial")
  if (ncol(coordinates(xy)) > 2) 
    stop("xy should be defined in two dimensions")
  as.what <- match.arg(as)
  prj <- proj4string(xy)
  ext <- extent(xy)
  # widen a tad (20%) the extent
  xplus <- (ext@xmax-ext@xmin) * 1.1
  yplus <- (ext@ymax-ext@ymin) * 1.1
  extplus <- extent(c(ext@xmin-xplus,ext@xmax+xplus,ext@ymin-yplus,ext@ymax+yplus)) 
  xsize <- ceiling((extplus@xmax-extplus@xmin) / res)
  ysize <- ceiling((extplus@ymax-extplus@ymin) / res)
  grid <- raster(extplus, nrows=ysize, ncols=xsize, crs=CRS(prj))
  res(grid) <- res
  values(grid) <- 1
  if(as.what=="SpatialPixels") {
    grid_sp <- as(grid, "SpatialPixels")  
    invisible(grid_sp)
  } else {
    invisible(grid)
  }
}
