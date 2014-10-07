\name{neg_cluster_matrix10302}
\alias{neg_cluster_matrix10302}
\docType{data}
\title{
Knn classification result of the negative controls of plate 10302
}
\description{
This data set contains an integer data frame with two colums.
The first column cotains the wellnamuber of origin of each sample cell in the screen.
The second column contains the population classification for each object in the screen.}
\usage{data(neg_cluster_matrix10302)}
\format{
  The format is:
 int [1:80047, 1:2] 24 24 24 24 24 24 24 24 24 24 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "origin" "raw_levels"
}
\examples{
data(neg_cluster_matrix10302)
print(negclust)
}
\keyword{datasets}
