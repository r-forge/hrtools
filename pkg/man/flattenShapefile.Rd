\name{flattenShapefile}
\alias{flattenShapefile}
\title{Convert 3D Shapefiles into 2D}
\usage{
  flattenShapefile(dsn = getwd(), layer = NULL)
}
\arguments{
  \item{dsn}{data source name (the folder containing the 3D
  ESRI Shapefile)}

  \item{layer}{layer name (the Shapefile name without
  extension)}
}
\value{
  An ESRI Shapefile named as the input (\code{layer}), plus
  the "\code{2D}" suffix will be written in the same
  directory (\code{dsn}) of the source spahefile. A
  \code{SpatialPoint} or \code{SpatialPointDataFrame} is
  returned.
}
\description{
  \code{adehabitatHR} works with 2D point locations.
  Possibly, data in ESRI Shapefile format could contain
  elevation information, that causes trpuble. This function
  'flattens' 3D Shapefiles into a 2D version. A 2D ESRI
  Shapefile will be written in the same directory
  (\code{dsn}) containing the input shapefile. The output
  shapefile will be named as the input (\code{layer}), plus
  the "\code{2D}" suffix.
}
\examples{
flat <- flattenShapefile(dsn='~/shapes/', layer='squirrels')
}
\references{
  \link{http://www.gdal.org/ogr/},
  \link{http://www.gdal.org/ogr/ogr_formats.html}
}
\seealso{
  \code{\link{readOGR}} and \code{\link{writeOGR}}
}

