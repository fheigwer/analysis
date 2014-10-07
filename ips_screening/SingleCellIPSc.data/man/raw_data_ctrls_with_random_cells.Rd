\name{raw_data_ctrls_with_random_cells}
\alias{raw_data_ctrls_with_random_cells}
\docType{data}
\title{
Raw control data set containing all features un filtered
}
\description{
This dataset contains the entire control object feature vectors aqcuired by CellProfiler and read in R.
All objects from each control well in the screen (positive and negative controls) were taken alongside with every 50 th other cell of the screen.
}
\usage{data(raw_data_ctrls_with_random_cells)}
\format{
  A data frame with 919796 observations on the following 67 variables.
  \describe{
    \item{\code{"ImageNumber"}}{a numeric vector}
\item{\code{"ObjectNumber"}}{a numeric vector}
\item{\code{"AreaShape_Area"}}{a numeric vector}
\item{\code{"AreaShape_Center_X"}}{a numeric vector}
\item{\code{"AreaShape_Center_Y"}}{a numeric vector}
\item{\code{"AreaShape_Compactness"}}{a numeric vector}
\item{\code{"AreaShape_Eccentricity"}}{a numeric vector}
\item{\code{"AreaShape_EulerNumber"}}{a numeric vector}
\item{\code{"AreaShape_Extent"}}{a numeric vector}
\item{\code{"AreaShape_FormFactor"}}{a numeric vector}
\item{\code{"AreaShape_MajorAxisLength"}}{a numeric vector}
\item{\code{"AreaShape_MinorAxisLength"}}{a numeric vector}
\item{\code{"AreaShape_Orientation"}}{a numeric vector}
\item{\code{"AreaShape_Perimeter"}}{a numeric vector}
\item{\code{"AreaShape_SolidityIntensity_IntegratedIntensityEdge_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_IntegratedIntensityEdge_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_IntegratedIntensityEdge_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_IntegratedIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_IntegratedIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_IntegratedIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_LowerQuartileIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_LowerQuartileIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_LowerQuartileIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MassDisplacement_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MassDisplacement_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MassDisplacement_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensityEdge_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensityEdge_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensityEdge_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MaxIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MeanIntensityEdge_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MeanIntensityEdge_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MeanIntensityEdge_OrigReIntensity_MeanIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MeanIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MeanIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MedianIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MedianIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MedianIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_MinIntensityEdge_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MinIntensityEdge_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MinIntensityEdge_OrigRedIntensity_MinIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_MinIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_MinIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_StdIntensityEdge_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_StdIntensityEdge_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_StdIntensityEdge_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_StdIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_StdIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_StdIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Intensity_UpperQuartileIntensity_OrigBlue"}}{a numeric vector}
\item{\code{"Intensity_UpperQuartileIntensity_OrigGreen"}}{a numeric vector}
\item{\code{"Intensity_UpperQuartileIntensity_OrigRed"}}{a numeric vector}
\item{\code{"Location_Center_X"}}{a numeric vector}
\item{\code{"Location_Center_Y"}}{a numeric vector}
\item{\code{"Neighbors_AngleBetweenNeighbors_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_FirstClosestDistance_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_FirstClosestObjectNumber_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_NumberOfNeighbors_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_PercentTouching_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_SecondClosestDistance_Expanded"}}{a numeric vector}
\item{\code{"Neighbors_SecondClosestObjectNumber_Expanded"}}{a numeric vector}
\item{\code{"Number_Object_Number"}}{a numeric vector}
  }
}
\examples{
data(raw_data_ctrls_with_random_cells)
}
\keyword{datasets}
