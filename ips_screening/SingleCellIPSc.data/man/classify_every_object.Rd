\name{classify_every_object}
\alias{classify_every_object}
\title{
K nearest neighbour classify each non-control cell of the screen
}
\description{
This function serevs as a rapper for the k nearest neighbour classification of the cells of plate 10301 as an example.
It uses the results of the k-means clustering and its underlying data as training dataset for a k nearest neighbor classification.  
}
\usage{
classify_every_object(plate = "10301")
}
\arguments{
  \item{plate}{
  non needed
  }
}
\value{
returns and saves a object of class knn
}

\author{
Florian Heigwer
}
\note{
k=100
}

\seealso{
\code{\link{knn}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (plate = "10301") 
{
    data("kmeans_result", package = "SingleCellIPSc.data")
    data("raw_controls", package = "SingleCellIPSc.data")
    data("raw_data_sample_10301", package = "SingleCellIPSc.data")
    data("origin10301", package = "SingleCellIPSc.data")
    clustering = cbind(origin, knn(raw_data, sample_data, raw_clust$cluster, 
        k = 100, prob = F))
    save(clustering, file = paste(plate, "clustering.RData", 
        sep = "_"))
  }
}
