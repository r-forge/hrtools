\name{squirrels}
\alias{squirrels}
\docType{data}
\title{
Radio-tracking Eurasian Red Squirrels
}
\description{
Data from Valtellina
}
\usage{data(squirrels)}
\format{
  The format is:
Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  ..@ data       :'data.frame':	1133 obs. of  12 variables:
  .. ..$ ID__TAG_: num [1:1133] 2429 2429 2429 2429 2429 ...
  .. ..$ AGE     : num [1:1133] 3 3 3 3 3 3 3 3 3 3 ...
  .. ..$ SEX     : num [1:1133] 2 2 2 2 2 2 2 2 2 2 ...
  .. ..$ DAY     : num [1:1133] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
  .. ..$ MO      : num [1:1133] 6 6 6 6 6 6 6 6 6 7 ...
  .. ..$ YEAR    : num [1:1133] 2006 2006 2006 2006 2006 ...
  .. ..$ ACTIV   : num [1:1133] 1 2 2 2 2 1 2 2 2 4 ...
  .. ..$ TREE    : num [1:1133] 10 10 10 10 10 10 10 10 10 10 ...
  .. ..$ HOUR    : Factor w/ 454 levels "08.40","09.00",..: 265 85 220 55 331 124 283 287 296 92 ...
  .. ..$ ALTITUDE: num [1:1133] 1978 1969 1955 1964 1970 ...
  .. ..$ METEO   : num [1:1133] 1 1 1 1 1 2 1 1 2 1 ...
  .. ..$ WIND    : num [1:1133] 1 2 1 1 1 1 2 2 2 1 ...
  ..@ coords.nrs : num(0) 
  ..@ coords     : num [1:1133, 1:2] 1596097 1596169 1595993 1596028 1596168 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  ..@ bbox       : num [1:2, 1:2] 1595282 5155574 1596900 5157519
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. ..$ : chr [1:2] "min" "max"
  ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slots
  .. .. ..@ projargs: chr "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +units=m +no_defs"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(squirrels)
## maybe str(squirrels) ; plot(squirrels) ...
}
\keyword{datasets}
