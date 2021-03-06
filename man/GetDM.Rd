% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetDM.R
\name{GetDM}
\alias{GetDM}
\title{GetDM}
\usage{
GetDM(pao, cov_list, is.het, is.annual, psi = FALSE, coord.p = TRUE,
  coord.th = TRUE)
}
\arguments{
\item{pao}{.pao file for the species of interest}

\item{cov_list}{List containing the index of covariates to include in each model}

\item{is.het}{Logical; Should heterogeniety in detection be modeled? (default = TRUE)}

\item{is.annual}{Logical; Should detection vary annually? (default = TRUE)}

\item{coord.p}{Logical; Should detection vary by lat/long (default = TRUE)}

\item{coord.th}{Logical; Should local availabilty vary by lat/long (default = TRUE)}

\item{cov_names}{Vector containing the names of all covariates}

\item{psi_list}{List containing the index of covariates on psi to include in each model}

\item{index}{Integer indicating which element in cov_list/psi_list to use}
}
\value{
A list containing 6 design matrices

dm1 Design matrix for psi, theta, theta1

dm2 Design matrix for gamma

dm3 Design matrix for epsilon

dm4 Design matrix for p

dm5 Design matrix for theta0

dm6 Design matrix for p.pi
}
\description{
Create design matrices for Presence
}

