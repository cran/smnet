\name{simulateNetwork}
\alias{simulateNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate river network topology and data
}
\description{
Simulates a graph, where nodes correspond to stream segments and an edge denotes that two stream segments are connected via a confluence.  In addition data is simulated on the network, with one datum associated with each stream segment.  Not of practical use, but good for checking that other utilities in SmoothNetwork work.  Uses the igraph package to construct the graph.
}
\usage{
simulateNetwork(n.segments, beta = NULL, lambda = 1, subs = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.segments}{
Number of nodes the simulated network should have
}
  \item{lambda}{
Positive real number: level of smoothness the simulated data should have.  Small values mean less smoothness.  Default = 1.
}
\item{beta}{
Numeric vector of covariate effects.  If supplied, independent random normal covariate data is generated on each stream segment with true coefficients as in \code{beta}.
}
  \item{subs}{Number between 0 and 1.  Determines whether data is simulated on a random subset of the network nodes of size \code{subs*n.segments}.
}
}

\value{
\item{response}{Vector of simulated data of length \code{n.segments}, or if \code{subs} is supplied, of length \code{subs*n.segments}}
\item{response.locs}{Vector indicating which node each value of \code{response}came from}
\item{alpha}{Random intercept term, generated from U(0,1) draw}
\item{z}{Vector of true spatial effects (\code{response} - \code{alpha} - random noise)}
\item{Z}{If \code{subs} is supplied these are the true values for all nodes}
\item{wgts}{Vector of (flow) weights, based on Shreve order}
\item{adjacency}{\code{n.segments}*\code{n.segments} adjacency matrix of class \code{spam}}
\item{coords}{\code{n.segments}*2 matrix of coordinates of node locations, useful for plotting}
}

\references{
Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.sf.net
}
\author{
Alastair Rushworth
}

\seealso{
\code{\link{plot.smnet}}, \code{\link{summary.smnet}}
}
\examples{

# Generate some simulated data with 100 stream segments
x<-simulateNetwork(n.segments = 100, lambda = 1)

lmSim<-smnet(formula=response~1 + 
               network(adjacency = x$adjacency, weight = x$wgts,
               locs = x$response.locs), 
               data.object = x$obsdata,
               method = "AICC")

plot(lmSim, type = "nodes", coords = x$obsdata[,2:3])

}

