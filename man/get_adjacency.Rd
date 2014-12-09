\name{get_adjacency}
\alias{get_adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct an Adjacency Matrix}
\description{Builds a sparse adjacency matrix from a user specified \code{SSN} data directory, by extracting and processing the binaryID.db table.  The resulting output of this function is typically a required input for fitting additive network models to \code{SSN} objects using the main \link{smnet} function.}
\usage{
get_adjacency(ssn_directory, net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ssn_directory}{Required character string indicating the path to the location of the .ssn directory which contains the binaryID.db table}
   \item{net}{Integer specifying the particular stream network of interest within the \code{SSN} object.  Defaults to 1.}
}

\value{Square sparse matrix of class "\code{spam}" with row and column dimension equal to the number of stream segments.  If the i^th column has non-zero elements \eqn{j_1}{j_1} and \eqn{j_2}{j_2} then this indicates that \eqn{j_1}{j_1} and \eqn{j_2}{j_2} are direct upstream neighbours of i.  

If the \eqn{i^\textrm{th}}{i^th} column has sum 1, then this indicates that \eqn{i}{i} has only one upstream neighbour, and therefore no confluence lies between them; by default the spatial penalties treat these differently.
}

\author{Alastair Rushworth}



\seealso{
\code{\link{smnet}}
}
\examples{
# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
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





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{network}
\keyword{sparse}
\keyword{P-spline}