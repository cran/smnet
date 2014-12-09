\name{predict.smnet}
\alias{predict.smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict From a Stream Network Model.}
\description{
Get predictions and standard errors at a fixed set of spatial locations and covariate values from a model fitted by \code{\link{smnet}}.
}
\usage{\method{predict}{smnet}(object, newdata =  NULL,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{smnet}, usually the result of a call to smnet.
}
\item{newdata}{New design matrix at which to make predictions, not required if the input data to \code{\link{smnet}} was an \code{SSN} object with an associated set of prediction points.
}
\item{...}{other arguments passed to predict.smnet
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{predictions}{Vector of predictions corresponding to prediction points in the \code{SSN} input object}
\item{predictions.se}{Vector of prediction standard errors}
%% ...
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

# fit the model using smnet
lmP<-smnet(formula = Sim_Values~1 + 
          network(adjacency = adjacency, weight = "addfunccol", netID = 1), 
          data.object = sim1.ssn)
                    
plot(lmP, type = "segments")
plot(lmP, type = "nodes")

# make some predictions
z<-predict(lmP)
}

