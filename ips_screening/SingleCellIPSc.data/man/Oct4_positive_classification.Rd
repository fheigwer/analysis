\name{Oct4_positive_classfication}
\alias{Oct4_positive_classfication}
\docType{data}
\title{
classification of each positive control
}
\description{
This data set contains an integer data frame with two colums.
The first column cotains the wellnamuber of origin of each positive control cell in the screen.
The second column contains the population classification for each positive control object in the screen.
}
\usage{data(Oct4_positive_classfication)}
\format{
  The format is:
 int [1:125190, 1:2] 23 23 23 23 23 23 23 23 23 23 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "origin" "raw_levels"
}
\examples{
data(Oct4_positive_classfication)
print(posclust)
}
\keyword{datasets}
