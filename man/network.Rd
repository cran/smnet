\name{network}
\alias{network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specify Network Smoother in Formulae}
\description{
This function specifies all of the information required to smooth parameters over the segments of a stream network using an adjacency matrix, and a vector of flow weights.  
}
\usage{
network(adjacency = NULL, weight = NULL, netID = 1, locs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{adjacency}{A sparse adjacency matrix describing the flow connections across the network.  \code{adjacency} must be of class "\code{spam}", and is typically obtained from a call to \link{get_adjacency}.
}
\item{weight}{Either a numeric vector of flow weights or a character string indicating the column name of a numeric vector of flow weights contained in the \code{data.object} that has been passed to \link{smnet}.  When a vector of flow weights is provided, this must have length equal to the number of rows of \code{adjacency} and respect the same ordering.
}
\item{netID}{Integer identifying the specific network of interest within the \code{SSN} object supplied to \link{smnet}.  Ignored if the input data to \code{smnet} is a \code{data.frame} and not an \code{SSN} object.
}
\item{locs}{Either a character string referring to a column of stream segment locations inside the data.object that has been passed to \code{smnet}.  Otherwise a numeric vector of stream segment locations, these should be arranged so that they correspond to the row and column ordering of the adjacency matrix.  Ignored if the input data to \code{smnet} is a \code{data.frame} and not an \code{SSN} object.
}
}

\value{
A list combining the processed input components above.  For internal use within \code{smnet}.
\item{adjacency}{Sparse adjacency matrix}
\item{weight}{Numeric vector of flow weights}
\item{netID}{Integer identifying network of interest}
\item{locs}{Vector of stream segments locations}
}

\author{Alastair Rushworth}

\seealso{\link{smnet}, \link{get_adjacency}}
\examples{
# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
library(SSN)
set.seed(101)
## simulate a SpatialStreamNetwork object
raw1.ssn <- createSSN(n = 100,
    obsDesign = binomialDesign(50), predDesign = binomialDesign(50),
  importToR = TRUE, path = paste(tempdir(),"/sim1", sep = ""))

## create distance matrices, including between predicted and observed
createDistMat(raw1.ssn, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
raw1DFobs <- getSSNdata.frame(raw1.ssn, "Obs")
raw1DFpred <- getSSNdata.frame(raw1.ssn, "preds")

## add a continuous covariate randomly
raw1DFobs[,"X1"] <- rnorm(length(raw1DFobs[,1]))
raw1DFpred[,"X1"] <- rnorm(length(raw1DFpred[,1]))

## add a categorical covariate randomly
raw1DFobs[,"F1"] <- as.factor(sample.int(3,length(raw1DFobs[,1]), replace = TRUE))
raw1DFpred[,"F1"] <- as.factor(sample.int(3,length(raw1DFpred[,1]), replace = TRUE))

## simulate Gaussian data
sim1.out <- SimulateOnSSN(raw1.ssn,
  ObsSimDF = raw1DFobs,
	PredSimDF = raw1DFpred,
	PredID = "preds",
	formula = ~ X1 + F1,
	coefficients = c(1, .5, -1, 1),
	CorModels = c("Exponential.tailup", "Exponential.taildown"),
	use.nugget = TRUE,
	use.anisotropy = FALSE,
	CorParms = c(2, 5, 2, 5, 0.1),
	addfunccol = "addfunccol")

## extract the ssn.object
sim1.ssn <- sim1.out$ssn.object

## extract the observed and predicted data frames, now with simulated values
sim1DFobs <- getSSNdata.frame(sim1.ssn, "Obs")
sim1DFpred <- getSSNdata.frame(sim1.ssn, "preds")

## store simulated prediction values, and then create NAs in their place
sim1preds <- sim1DFpred[,"Sim_Values"]
sim1DFpred[,"Sim_Values"] <- NA
sim1.ssn <- putSSNdata.frame(sim1DFpred, sim1.ssn, "preds")

# create the adjacency matrix for use with SmoothNetwork
adjacency<-get_adjacency(paste(tempdir(),"/sim1", sep = ""))

lmP<-smnet(formula = Sim_Values~1 + 
          network(adjacency = adjacency, weight = "addfunccol", netID = 1), 
          data.object = sim1.ssn)
          
plot(lmP, type = "segments")
plot(lmP, type = "nodes")
}
