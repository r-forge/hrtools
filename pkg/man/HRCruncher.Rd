\name{HRCruncher}
\alias{HRCruncher}
\title{Batch home range calculations.}
\usage{
  HRCruncher(xy, method = c("mcp", "href", "lscv", "hadj"),
    grid = 20, same4all = FALSE, sizeGrid = TRUE,
    idfield = c("ID"), minfix = 15, percent = 95,
    keepfields = NULL, unin = c("m", "km"),
    unout = c("ha", "km2", "m2"), boundary = NULL, ...)
}
\arguments{
  \item{xy}{An object inheriting the class
  \code{SpatialPoints} containing the x and y relocations
  of the animal. If \code{xy} inherits the class
  \code{SpatialPointsDataFrame}, it should contain a column
  whose values specify the identity of the animals for each
  relocation. If that column exists and is not named
  \code{ID}, its name must be specified using the \code{id}
  argument. If a \code{SpatialPoints} object is passed, it
  is assumed that all the relocations are from a single
  animal.}

  \item{method}{a string specifying the home range
  calculation method to be used. See Details for further
  specifications.}

  \item{grid}{a number giving the \emph{cell size} (either
  as number of cells or in meters) of the grid on which the
  UD should be estimated. Alternatively, this parameter may
  be an object inheriting the class SpatialPixels, that can
  be created using the \code{sizeGrid} function.}

  \item{same4all}{logical. If \code{TRUE}, the same grid
  (see above) is used for all animals, if \code{FALSE}, for
  each animal a reference grid will be calculated.}

  \item{sizeGrid}{if \code{TRUE} (default) the UD
  estimation grid (see \code{grid} and \code{same4all}
  above) will be calculated using the \code{sizeGrid}
  function, instead of partitioning the extent in a
  \code{grid} x \code{grid} cells mesh, as "regular"
  \code{adehabitatHR} kernels do.}

  \item{idfield}{the name of a column in \code{xy}
  corresponding to the identity of the animals for each
  relocation.}

  \item{minfix}{an integer specifying the minimum number of
  relocations needed to calculate an home range. Default
  value is 15. If an animal has less than \code{minfix}
  relocations a warning will be issued and no home renge
  calculations will be made for that animal.}

  \item{percent}{for \code{mcp} method, the number of
  relocations to use to calculate a Minimum Convex Polygon.
  The default value of 95 means that a 5% of points falling
  far from the centroid of the relocations used will be
  treated as outliers. See \code{percent} argument of
  \code{mcp} function. For other (kernel-based) methods
  \code{prc} represents the percentage level (default: 95%)
  for home range estimation (see \code{percent} parameter
  in \code{kernelUD} function).}

  \item{keepfields}{a character vector (as from
  \code{names}) with names of fields to be preserved in the
  results. \code{idfield} is automatically kept. In case
  some field has different values, the first appearing
  value will be kept.}

  \item{unin}{the units of the relocations coordinates.
  Either "m" for meters (default) or "km" for kilometers.}

  \item{unout}{the units of the output areas. Either "m2"
  for square meters, "km2" for square kilometers or "ha"
  for hectares (default).}

  \item{boundary}{If not NULL, an object inheriting the
  class SpatialLines defining a barrier that cannot be
  crossed by the animals. See \code{kernelUD} for details.}

  \item{...}{other parameters, passed to \code{mcp} and/or
  \code{kernelUD}.}
}
\value{
  an object of \code{HRData} class, with its slots filled.
}
\description{
  Driver function for batch homerenge calculations.
}
\details{
  Available methods (see \code{method} above) are:
  \code{mcp}, \code{href}, \code{lscv}, \code{hadj}.
  \itemize{ \item{\code{mcp}}{ Minimum Convex Polygon,
  using \code{adehabitatHR} \code{mcp} function.}
  \item{\code{href}}{ kernel home range using
  \code{adehabitatHR} \code{kernelUD} function with "ad
  hoc" smoothing parameter \code{h}.} \item{\code{lscv}}{
  kernel home range using \code{adehabitatHR}
  \code{kernelUD} function with least-squares cross
  validation to estimate \code{h} smoothing parameter.}
  \item{\code{hadj}}{ kernel home range using
  \code{adehabitatHR} \code{kernelUD} function with
  "adjusted" \code{h} as in Wauters et al. (2007).} }
}
\examples{
data(squirrels)
# size up a grid with 30 m cell size
reference.grid <- sizeGrid(squirrels, res=30)

# plot reference grid and locations
image(reference.grid, col='lightgray')
plot(squirrels, col=seq(1:18), add=TRUE)

# define animal ID field
idfield <- "TAG"

# define minumum number of fixes needed
minfix <- 15

calculate mcp and hadj kernel (href and lscv kernels will be calculated implicitly)
result <- HRCruncher(squirrels, method=c('mcp','hadj'), percent=95, grid=reference.grid, idfield=idfield, unin='m', unout='ha')

# save results as shapefile
libaray(rgdal)
writeOGR(result@href@geometry, dsn='/tmp', layer='test' , driver="ESRI Shapefile")
}
\references{
  Lucas A. Wauters, Damiano G. Preatoni, Ambrogio Molinari,
  Guido Tosi (2007). Radio-tracking squirrels: Performance
  of home range density and linkage estimators with small
  range and sample size. Ecological Modelling 202(3-4):
  333-344
  \url{http://dx.doi.org/10.1016/j.ecolmodel.2006.11.001}
}
\seealso{
  \link{sizeGrid}, \link[adehabitatHR]{mcp},
  \link[adehabitatHR]{kernelUD}, \link{HRData-class}
}

