\name{sizeGrid}
\alias{sizeGrid}
\title{Create a reference grid for kernel calculations}
\usage{
  sizeGrid(xy, res = 100,
    as = c("SpatialPixels", "raster"))
}
\arguments{
  \item{xy}{an object inheriting the class
  \code{SpatialPoints} containing the x and y relocations
  of the animal. That object will be used to calculate grid
  extent and CRS.}

  \item{res}{The grid spatial resolution, in meters.}

  \item{as}{specify whether to output the reference grid as
  \code{SpatialPixels} object (default) or as
  \code{raster}.}
}
\value{
  a \code{SpatialPixels} object (default) or a
  \code{raster} object.
}
\description{
  Create a \code{SpatialPixels} object with same extent and
  CRS of input points and resolution as desired.
}
\examples{
# generate a reference grid with a 50 m spatial resolution for the squirrel dataset
data(squirrels)
reference.grid <- sizeGrid(squirrels, res=50)
}

