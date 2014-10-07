\name{clustering_sample_10303}
\alias{clustering_sample_10303}
\docType{data}
\title{
classification results for each object in plate 10303
}
\description{
The data object called clustering contains a data frame of two columns.
The clomuns contain the well of origin and the knn classification of each object in plate 10303.
}
\usage{data(clustering_sample_10303)}
\format{
  A data frame with x observations on the following 2 variables.
  \describe{
    \item{\code{1}}{a numeric vector}
    \item{\code{2}}{a numeric vector}
  }
}
\details{
the first column contains the origin of each object of plate 10303
the second column contains class of each object of plate 10303
}
\examples{
data(clustering_sample_10303)
hist(clustering[,2])
}
\keyword{datasets}
