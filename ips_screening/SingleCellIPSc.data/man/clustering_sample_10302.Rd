\name{clustering_sample_10302}
\alias{clustering_sample_10302}
\docType{data}
\title{
classification results for each object in plate 10302
}
\description{
The data object called clustering contains a data frame of two columns.
The clomuns contain the well of origin and the knn classification of each object in plate 10302.
}
\usage{data(clustering_sample_10302)}
\format{
  A data frame with x observations on the following 2 variables.
  \describe{
    \item{\code{1}}{a numeric vector}
    \item{\code{2}}{a numeric vector}
  }
}
\details{
the first column contains the origin of each object of plate 10302
the second column contains class of each object of plate 10302
}
\examples{
data(clustering_sample_10302)
hist(clustering[,2])
}
\keyword{datasets}
