\name{all_samples_annotation}
\alias{all_samples_annotation}
\docType{data}
\title{
Wellnames of all plates
}
\description{
This dataset contains the well count, well pace and RNAi reagent for each treatment on all plates (10301,10302,10303)
}
\usage{data(all_samples_annotation)}
\format{
  A data frame with 2112 observations on the following 3 variables.
  \describe{
    \item{\code{V1}}{a numeric vector}
    \item{\code{V2}}{a factor with levels \code{A03} \code{A04} \code{A05} \code{A06} \code{A07} \code{A08} \code{A09} \code{A10} \code{A11} \code{A12} \code{A13} \code{A14} \code{A15} \code{A16} \code{A17} \code{A18} \code{A19} \code{A20} \code{A21} \code{A22} \code{A23} \code{A24} \code{B03} \code{B04} \code{B05} \code{B06} \code{B07} \code{B08} \code{B09} \code{B10} \code{B11} \code{B12} \code{B13} \code{B14} \code{B15} \code{B16} \code{B17} \code{B18} \code{B19} \code{B20} \code{B21} \code{B22} \code{B23} \code{B24} \code{C03} \code{C04} \code{C05} \code{C06} \code{C07} \code{C08} \code{C09} \code{C10} \code{C11} \code{C12} \code{C13} \code{C14} \code{C15} \code{C16} \code{C17} \code{C18} \code{C19} \code{C20} \code{C21} \code{C22} \code{C23} \code{C24} \code{D03} \code{D04} \code{D05} \code{D06} \code{D07} \code{D08} \code{D09} \code{D10} \code{D11} \code{D12} \code{D13} \code{D14} \code{D15} \code{D16} \code{D17} \code{D18} \code{D19} \code{D20} \code{D21} \code{D22} \code{D23} \code{D24} \code{E03} \code{E04} \code{E05} \code{E06} \code{E07} \code{E08} \code{E09} \code{E10} \code{E11} \code{E12} \code{E13} \code{E14} \code{E15} \code{E16} \code{E17} \code{E18} \code{E19} \code{E20} \code{E21} \code{E22} \code{E23} \code{E24} \code{F03} \code{F04} \code{F05} \code{F06} \code{F07} \code{F08} \code{F09} \code{F10} \code{F11} \code{F12} \code{F13} \code{F14} \code{F15} \code{F16} \code{F17} \code{F18} \code{F19} \code{F20} \code{F21} \code{F22} \code{F23} \code{F24} \code{G03} \code{G04} \code{G05} \code{G06} \code{G07} \code{G08} \code{G09} \code{G10} \code{G11} \code{G12} \code{G13} \code{G14} \code{G15} \code{G16} \code{G17} \code{G18} \code{G19} \code{G20} \code{G21} \code{G22} \code{G23} \code{G24} \code{H03} \code{H04} \code{H05} \code{H06} \code{H07} \code{H08} \code{H09} \code{H10} \code{H11} \code{H12} \code{H13} \code{H14} \code{H15} \code{H16} \code{H17} \code{H18} \code{H19} \code{H20} \code{H21} \code{H22} \code{H23} \code{H24} \code{I03} \code{I04} \code{I05} \code{I06} \code{I07} \code{I08} \code{I09} \code{I10} \code{I11} \code{I12} \code{I13} \code{I14} \code{I15} \code{I16} \code{I17} \code{I18} \code{I19} \code{I20} \code{I21} \code{I22} \code{I23} \code{I24} \code{J03} \code{J04} \code{J05} \code{J06} \code{J07} \code{J08} \code{J09} \code{J10} \code{J11} \code{J12} \code{J13} \code{J14} \code{J15} \code{J16} \code{J17} \code{J18} \code{J19} \code{J20} \code{J21} \code{J22} \code{J23} \code{J24} \code{K03} \code{K04} \code{K05} \code{K06} \code{K07} \code{K08} \code{K09} \code{K10} \code{K11} \code{K12} \code{K13} \code{K14} \code{K15} \code{K16} \code{K17} \code{K18} \code{K19} \code{K20} \code{K21} \code{K22} \code{K23} \code{K24} \code{L03} \code{L04} \code{L05} \code{L06} \code{L07} \code{L08} \code{L09} \code{L10} \code{L11} \code{L12} \code{L13} \code{L14} \code{L15} \code{L16} \code{L17} \code{L18} \code{L19} \code{L20} \code{L21} \code{L22} \code{L23} \code{L24} \code{M03} \code{M04} \code{M05} \code{M06} \code{M07} \code{M08} \code{M09} \code{M10} \code{M11} \code{M12} \code{M13} \code{M14} \code{M15} \code{M16} \code{M17} \code{M18} \code{M19} \code{M20} \code{M21} \code{M22} \code{M23} \code{M24} \code{N03} \code{N04} \code{N05} \code{N06} \code{N07} \code{N08} \code{N09} \code{N10} \code{N11} \code{N12} \code{N13} \code{N14} \code{N15} \code{N16} \code{N17} \code{N18} \code{N19} \code{N20} \code{N21} \code{N22} \code{N23} \code{N24} \code{O03} \code{O04} \code{O05} \code{O06} \code{O07} \code{O08} \code{O09} \code{O10} \code{O11} \code{O12} \code{O13} \code{O14} \code{O15} \code{O16} \code{O17} \code{O18} \code{O19} \code{O20} \code{O21} \code{O22} \code{O23} \code{O24} \code{P03} \code{P04} \code{P05} \code{P06} \code{P07} \code{P08} \code{P09} \code{P10} \code{P11} \code{P12} \code{P13} \code{P14} \code{P15} \code{P16} \code{P17} \code{P18} \code{P19} \code{P20} \code{P21} \code{P22} \code{P23} \code{P24}}
    \item{\code{V3}}{a character vector}
  }
}
\details{
the firt column contains the well numbers see(origin)
the second column contains the well names see(origin)
the third column contains the RNAi eagent used in that well
}
\examples{
data(all_samples_annotation)
print(names_complete)
}
\keyword{datasets}
