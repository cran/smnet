% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.smnet.R
\name{predict.smnet}
\alias{predict.smnet}
\title{Make Predictions Using a Fitted Stream Network Model.}
\usage{
\method{predict}{smnet}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{Object of class \code{smnet}, usually the result of a call to \code{smnet()}}

\item{newdata}{\code{data.frame} containing new covariate values at which to make predictions.}

\item{...}{other arguments passed to \code{predict.smnet()}}
}
\value{
List object with components
\itemize{
\item{\code{predictions}: vector of predictions corresponding to prediction points in
the original \code{SpatialStreamNetwork} input object}
\item{\code{predictions.se}: vector of prediction standard errors}
}
}
\description{
Get predictions and standard errors for a new set of spatial locations
and associated covariate values from a model fitted by \code{smnet()}.
}
\author{
Alastair Rushworth
}
