\name{plot.smnet}
\alias{plot.smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a Stream Network Model}
\description{Plot linear, univariate and bivariate smooth effects and network smooth terms that resulting from a call to \code{\link{smnet}}.  }
\usage{
  \method{plot}{smnet}(x, type = "covariates", se = FALSE, res = FALSE, coords = NULL, key = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \bold{smnet}
  }
  \item{type}{Character string identifying the type of plot to be produced.  The default, "\code{covariates}", produces plots of all linear and smooth components (the latter corresponding to each appearance of \code{\link{m}} in the model formula).   "\code{node}" plots the spatial network model component associated with a use of \code{\link{network}} in the model formula; the plot uses a set of nodes each associated with a stream stretch and plotted at their geographical midpoints.  The colour of each node corresponds to the fitted value associated with that stream segment.  "\code{segments}" plots the spatial network model component associated with a use of \code{\link{network}} in the model formula; segments are plotted using information within an \bold{SSN} object (may be slower than "\code{node}" for larger networks).
  }
  \item{se}{Logical.  When TRUE (the default), calculates and adds standard errors to plots of linear and smooth components the default is TRUE.  Setting to FALSE may be be quicker for very large data sets.}
  \item{res}{Logical.  Plots partial residuals on linear and smooth component plots when TRUE, the default.  Ignored if cov = FALSE.}
  \item{coords}{A 2-column matrix of coordinates required to plot spatial smooths when input data was not an \bold{SSN} object
  }
  \item{key}{Logical.  Plots a colour legend for the node plot when set to TRUE (the default.  Ignored when node = FALSE.)
  }
    \item{...}{Other arguments passed to plot
  }
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

lmP<-smnet(formula = Sim_Values~1 + 
          network(adjacency = adjacency, weight = "addfunccol", netID = 1), 
          data.object = sim1.ssn)
          
plot(lmP, type = "segments")
plot(lmP, type = "nodes")
}


\author{Alastair Rushworth}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
   \code{\link{predict.smnet}}, \code{\link{summary.smnet}}
  }


