\name{make_and_save_kmeans}
\alias{make_and_save_kmeans}
\title{
K-means clsuter the control data set
}
\description{
This function is used to cluster the crontrol data set into 6 distinct phenotypic groups and save results accordingly.

30 random starts were used in this analysis with k=6.
}
\usage{
make_and_save_kmeans()
}
\arguments{
  NaN
}
\details{
random starts=30
k=6
}
\value{
saves the clsutering result and the data matrix as the RData objects raw_clust and raw_data
}
\author{
Florian Heigwer
}
\note{
This functions takes a very long time (h)
}
\seealso{
\code{\link{kmeans}}
}
\examples{
##do not run, it will take a very long time 

#make_and_save_kmeans()

## The function is currently defined as
function () 
{
    data("raw_data_ctrls_with_random_cells",package="SingleCellIPSc.data")
    drop_v <- c(1, 2, 4, 5, 58, 59, 67, 62, 66, 16, 17, 18, 25, 
        26, 27, 8, 10)
    keep <- c(1, 4, 20, 19, 25, 26, 47, 48)
    raw_data <- raw_data_ctrls_with_random_cells[, c(keep)]
    for (y in 1:dim(raw_data)[2]) {
        raw_data[, y] = (raw_data[, y] - min(raw_data[, y]))/(max(raw_data[, 
            y]) - min(raw_data[, y]))
    }
    raw_clust <- kmeans(raw_data, 6, nstart = 30)
    save(raw_data, file = "raw_controls.RData")
    save(raw_clust, file = "kmeans_result.RData")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
