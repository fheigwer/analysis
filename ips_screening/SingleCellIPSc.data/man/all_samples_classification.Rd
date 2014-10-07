\name{all_samples_classifications}
\alias{all_samples_classifications}
\docType{data}
\title{
Knn classification result of all sample cells
}
\description{
This data set contains an integer data frame with two colums.
The first column cotains the wellnamuber of origin of each sample cell in the screen.
The second column contains the population classification for each object in the screen.
}
\usage{data(all_samples_classifications)}
\format{
  The format is:
 int [1:7731725, 1:2] 3 3 3 3 3 3 3 3 3 3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "origin" "raw_levels"
}
\examples{
data(all_samples_classifications)
print(clustering_complete)
}
\keyword{datasets}
