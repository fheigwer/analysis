\name{SingleCellIPSc.data-package}
\alias{SingleCellIPSc.data-package}
\alias{SingleCellIPSc.data}
\docType{package}
\title{
This package bundles all data and scripts used in the preparation of the manuscript:

Volz et al., Single-cell phenotyping of human iPSc using image-based 
high throughput screening

}
\description{
We adapted conditions that allow for the utilization of hiPSc in RNA interference (RNAi) screening and enhanced image analysis by single cell phenotyping to account for the inherently increased phenotypic variability of these cells. Using our newly established high-throughput RNAi assay and single cell feature based classification we were able to identify, characterize and validate key regulators of pluripotency in hiPSc formerly masked in “aggregate features per well analysis”.

Therefor we provide this package with simple functions and wrappers for everybody, who wants to use this method on his own data sets.
}
\details{
\tabular{ll}{
Package: \tab SingleCellIPSc.data\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2014-07-22\cr
License: \tab GPL-2\cr
}

}
\author{
Florian Heigwer, Christian Volz, Michael Boutros

Maintainer: Who to complain to <f.heigwer@dkfz.de>
}
\references{
Volz et al. NMEHT under submission
}

\keyword{ package }
\examples{
make_scatter_plot_samples(sample_name = "TAF1", size = 2000, x = 7, y = 1, z = 6, alpha = 0.4)
make_scatter_plot_samples(size = 2000, x = 7, y = 1, z = 6, alpha = 0.4)
## make_distance_plot(absolute = FALSE, relative = TRUE, feature = FALSE)
}
