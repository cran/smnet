\name{smnet}
\alias{smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Additive Modelling for Stream Networks}

\description{Fits (Gaussian) additive models to river network data based on the flexible modelling framework described in O'Donnell et al. (2014).  Data must either be in the form of an object of class \code{SpatialStreamNetwork} as used by the R package \code{SSN}  (Ver Hoef et al., 2012) or in the form of a \code{data.frame}.  Smoothness of covariate effects is represented by transforming data onto a uniformly spaced B-spline basis; parameter estimates are obtained by penalized least-squares.  Optimal smoothness is achieved by using a numerical optimization of \code{AIC}, \code{GCV} or \code{AICC}.  

Care is taken to exploit the sparse matrix properties of model objects using the optimized storage and algebra routines in the R packages \code{spam} (Furrer and Sain, 2010).   This allows faster fitting and lower memory footprint.

The formula interpreter used for penalised additive components is modelled on the code found in the package \code{mgcv}. 
}

\usage{
smnet(formula, data.object, control = NULL, method = "AICC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula similar to those as used in by the gam function in the package \code{mgcv}.  Smooth functions based on P-splines with \code{m(..., k=20)} function, up to 2-dimensional interactions currently supported.  At present, only \code{k} the spline basis size can be specified using \link{m}.}
  \item{data.object}{Either an object of class "SpatialStreamNetwork" or a \code{data.frame} containing response variable and covariates}
\item{control}{A list of options passed to the optimiser.  \code{maxit}, default = 500, sets an upper limit of iterations made by the optimiser.  \code{approx  = NULL}, positive integer specifying the number samples to collect using a Monte-Carlo method to approximate the performance criterion surface when direct evaluation is too slow - this takes advantage of matrix sparsity and may be much faster if the network has a large number of segments or the data is large (\code{approx  = 100} often works well in these cases.}
  \item{method}{Character string determining the performance criterion for choosing optimal smoothness, options are \code{"AICC"} or \code{"GCV"}. 
}
}

\value{Object of class \code{smnet} with components
\item{1}{Original \code{SSN} object used provided as \code{data.object}, unchanged}
\item{2}{List; model output including fitted values (\code{fit}), effective dimension (\code{ED}), residual variance (\code{sigma.sq}),  sparse matrices involved in model fit etc.
}}
\references{
Ver Hoef, J.M.,  Peterson, E.E., Clifford, D., Shah, R.  (2012)   SSN: An R Package for Spatial Statistical Modeling on Stream Networks 

O' Donnell, D., Rushworth, A.M., Bowman, A.W., Scott, E.M., Hallard, M.  (2014) Flexible regression models over river networks.  Journal of the Royal Statistical Society: Series C (Applied Statistics). 63(1) 47--63.

Reinhard Furrer, Stephan R. Sain (2010). spam: A Sparse Matrix R
Package with Emphasis on MCMC Methods for Gaussian Markov Random
Fields. Journal of Statistical Software, 36(10), 1-25. URL: http://www.jstatsoft.org/v36/i10/
}
\author{Alastair Rushworth}


\seealso{
\code{\link[=smnet]{get_adjacency}}, \code{\link{plot.smnet}}
}

\examples{
# Set up an SSN object; this part is taken from the SSN:::SimulateOnSSN help file
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

lmP<-smnet(formula = Sim_Values~1 + m(NEAR_X)+
          network(adjacency = adjacency, weight = "addfunccol", netID = 1), 
          data.object = sim1.ssn)
          
plot(lmP, type = "covariates")
plot(lmP, type = "nodes")
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{network}
\keyword{sparse}
\keyword{P-spline}