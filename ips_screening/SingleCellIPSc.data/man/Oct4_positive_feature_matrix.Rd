\name{Oct4_positive_feature_matrix}
\alias{Oct4_positive_feature_matrix}
\docType{data}
\title{
positive control data for all plates
}
\description{
This data set contains all filtered feature data for the control objects of all plates.
The data was feature wise scaled from 0 to 1.
}
\usage{data(Oct4_positive_feature_matrix)}
\format{
  A data frame with 125190 observations on the following 8 variables.
  \describe{
    \item{\code{Nuclear Area}}{a numeric vector}
    \item{\code{Nuclear Extent}}{a numeric vector}
    \item{\code{max. Oct4 Intensity}}{a numeric vector}
    \item{\code{max. Hoechst Intensity}}{a numeric vector}
    \item{\code{ø Hoechst Intensity}}{a numeric vector}
    \item{\code{ø Oct4 Intensity}}{a numeric vector}
    \item{\code{Distance from first neighbor}}{a numeric vector}
    \item{\code{Number of Neighbors}}{a numeric vector}
  }
}
\examples{
data(Oct4_positive_feature_matrix)
print(posmatrix)
}
\keyword{datasets}
