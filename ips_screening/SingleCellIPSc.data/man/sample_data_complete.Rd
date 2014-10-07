\name{all_samples_feature_data}
\alias{all_samples_feature_data}
\docType{data}
\title{
combined sample data set containing all 8 filtered features 
}
\description{
This data set contains all filtered feature data for every object in the screen, which is not from the control treatments.
}
\usage{data(all_samples_feature_data)}
\format{
  A data frame with 7731725 observations on the following 8 variables.
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
\examples{
data(all_samples_feature_data)
pairs(sample_data_complete)
}
\keyword{datasets}
