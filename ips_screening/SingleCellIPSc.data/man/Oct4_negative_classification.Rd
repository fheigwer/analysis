\name{Oct4_negative_classfication}
\alias{Oct4_negative_classfication}
\docType{data}
\title{
classification of each negative control
}
\description{
This data set contains an integer data frame with two colums.
The first column cotains the wellnamuber of origin of each negative control cell in the screen.
The second column contains the population classification for each negative control object in the screen.
}
\usage{data(Oct4_negative_classfication)}
\format{
  The format is:
 int [1:234586, 1:2] 24 24 24 24 24 24 24 24 24 24 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "origin" "raw_levels"
}
\examples{
data(Oct4_negative_classfication)
print(negclust)
}
\keyword{datasets}
