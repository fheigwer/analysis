\name{make_scatter_plot_ctrls}
\alias{make_scatter_plot_ctrls}

\title{
Scatterplot cells of a sample in 3D
}
\description{
This function can be used to reproduce te scatterplots shown in figure 1c
}
\usage{
make_scatter_plot_ctrls(size = 20000, x = 7, y = 1, z = 6, alpha = 0.4)
}
\arguments{
  \item{size}{
    \code{integer}{how many dots should be plotted}
}
  \item{x}{
    \code{integer}{column to be the x axis}
}
  \item{y}{
 \code{integer}{column to be the y axis}
}
  \item{z}{
 \code{integer}{column to be the z axis}
}
  \item{alpha}{
    \code{float}{defines the opacity of points in the plot}
}
}
\value{
this funtion does not return a value but a plot
}
\references{
scatterplot3d
}
\author{
Florian Heigwer
}
\seealso{
\code{\link{make_scatter_plot_samples}},\code{\link{scatterplot3d}}
}
\examples{

make_scatter_plot_samples(size = 2000, x = 7, y = 1, z = 6, alpha = 0.4)

## The function is currently defined as
function (size = 20000, x = 7, y = 1, z = 6, 
    alpha = 0.4) 
{
    data("Oct4_positive_classfication", package = "SingleCellIPSc.data")
    data("Oct4_positive_feature_matrix", package = "SingleCellIPSc.data")
    data("Oct4_negative_classfication", package = "SingleCellIPSc.data")
    data("Oct4_negative_feature_matrix", package = "SingleCellIPSc.data")
    if (size > length(posclust[, 2])) {
        size <- length(posclust[, 2])
    }
    sample = sample(1:(length(posclust[, 2])), size = size)
    par(mfrow = c(1, 2))
    scatterplot3d(posmatrix[sample, x], posmatrix[sample, y], 
        posmatrix[sample, z], col.lab = "black", col.axis = "black", 
        col.grid = "white", box = FALSE, pch = ".", xlab = colnames(negmatrix)[7], 
        ylab = colnames(negmatrix)[1], zlab = colnames(negmatrix)[6], 
        xlim = c(-0.1, 0.5), ylim = c(0, 0.8), zlim = c(0, 0.6), 
        color = alpha(posclust[sample, 2], alpha))
    if (size > length(negclust[, 2])) {
        size <- length(negclust[, 2])
    }
    sample = sample(1:(length(negclust[, 2])), size = size)
    scatterplot3d(negmatrix[sample, x], negmatrix[sample, y], 
        negmatrix[sample, z], col.lab = "black", col.axis = "black", 
        col.grid = "white", box = FALSE, pch = ".", xlab = colnames(negmatrix)[7], 
        ylab = colnames(negmatrix)[1], zlab = colnames(negmatrix)[6], 
        xlim = c(-0.1, 0.5), ylim = c(0, 0.8), zlim = c(0, 0.6), 
        color = alpha(negclust[sample, 2], alpha))
  }
}
