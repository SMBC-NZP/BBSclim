% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetOcc.R
\name{GetOccProb}
\alias{GetOccProb}
\title{GetOccProb}
\usage{
GetOccProb(alpha, buff_method = "rec", Write = TRUE)
}
\arguments{
\item{alpha}{Four letter alpha code for species of interest}

\item{buff_method}{Optional: Which buffer method to use? Options include: "rec", "bcr"}

\item{buffer}{Dataframe containing the buffered spp count data}
}
\description{
Estimate annual probability of occupancy using parameter estimates from top model
}

