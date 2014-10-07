\name{kmeans_result}
\alias{kmeans_result}
\docType{data}
\title{
K means clustering result 
}
\description{
This data set contains the kmeans object given from the kmeans clustering of the raw_controls data set.
There are defined the 6 clsuter centers, which are the later population characterstics. Also this serves as the training data set for the knn classification of all other cells in the screen.
}
\usage{data(kmeans_result)}
\format{
  The format is:
List of 9
 $ cluster     : Named int [1:919644] 3 6 3 3 4 5 6 5 3 3 ...
  ..- attr(*, "names")= chr [1:919644] "1" "2" "3" "4" ...
 $ centers     : num [1:6, 1:8] 0.29 0.187 0.126 0.201 0.148 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:6] "1" "2" "3" "4" ...
  .. ..$ : chr [1:8] "AreaShape_Area" "AreaShape_Extent" "Intensity_MaxIntensity_OrigGreen" "Intensity_MaxIntensity_OrigBlue" ...
 $ totss       : num 143113
 $ withinss    : num [1:6] 9290 6018 7951 10014 8113 ...
 $ tot.withinss: num 47884
 $ betweenss   : num 95228
 $ size        : int [1:6] 93049 79850 239617 164948 214676 127504
 $ iter        : int 6
 $ ifault      : int 0
 - attr(*, "class")= chr "kmeans"
}
\examples{
data(kmeans_result)
print(raw_clust)
}
\keyword{datasets}
