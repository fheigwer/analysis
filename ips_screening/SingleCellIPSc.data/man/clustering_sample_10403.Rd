\name{clustering_sample_10403}
\alias{clustering_sample_10403}
\docType{data}
\title{
classification results for each object in plate 10403
}
\description{
The data object called clustering contains a data frame of two columns.
The clomuns contain the well of origin and the knn classification of each object in plate 10403.
}
\usage{data(clustering_sample_10403)}
\format{
  A data frame with x observations on the following 2 variables.
  \describe{
    \item{\code{1}}{a numeric vector}
    \item{\code{2}}{a numeric vector}
  }
}
\details{
the first column contains the origin of each object of plate 10403
the second column contains class of each object of plate 10403
}
\examples{
data(clustering_sample_10403)
hist(clustering[,2])
}
\keyword{datasets}
