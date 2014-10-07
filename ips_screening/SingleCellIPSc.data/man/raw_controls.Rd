\name{raw_controls}
\alias{raw_controls}
\docType{data}
\title{
combined control data set containing all 8 filtered features 
}
\description{
This data set contains all filtered feature data for every object in the screen, which is from the control treatments.
}
\usage{data(raw_controls)}
\format{
  A data frame with 919644 observations on the following 8 variables.
  \describe{
    \item{\code{AreaShape_Area}}{a numeric vector}
    \item{\code{AreaShape_Extent}}{a numeric vector}
    \item{\code{Intensity_MaxIntensity_OrigGreen}}{a numeric vector}
    \item{\code{Intensity_MaxIntensity_OrigBlue}}{a numeric vector}
    \item{\code{Intensity_MeanIntensity_OrigBlue}}{a numeric vector}
    \item{\code{Intensity_MeanIntensity_OrigGreen}}{a numeric vector}
    \item{\code{Neighbors_FirstClosestDistance_Expanded}}{a numeric vector}
    \item{\code{Neighbors_NumberOfNeighbors_Expanded}}{a numeric vector}
  }
}
\details{
This is a subset of the data from raw_data_ctrls_with_random_cells.
}
\examples{
data(raw_controls)
pairs(raw_data)
}
\keyword{datasets}
