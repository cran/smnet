\name{m}
\alias{m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify Smooth Terms in Formulae
}
\description{
Function used to set up univariate or bivariate smooth terms based on P-splines, for use within a call to \code{\link{smnet}}.   
}
\usage{
m(..., k = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
One or variables for which a P-spline smooth is desired
}
  \item{k}{
Integer, defines the number of evenly spaced B-spline basis functions to represent the smooth component, default is 10.  For two-dimensional smooths, this is the marginal basis size.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{term}{Character vector of the names of the variables involved in the smooth to be set up}
\item{bs.dim}{Number of B-spline basis functions to be used in the smooth}
}
\references{
Modified version of \code{s} originally from package \code{mgcv}, Simon Wood (2014).
}
\author{
Alastair Rushworth
}



\seealso{
\code{\link{smnet}}
}
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

lmP<-smnet(formula = Sim_Values~ m(NEAR_X) +
          network(adjacency = adjacency, weight = "addfunccol",
          netID = 1), data.object = sim1.ssn)
          
plot(lmP, type = "covariates", res = TRUE)      
plot(lmP, type = "segments")
plot(lmP, type = "nodes")      
}

