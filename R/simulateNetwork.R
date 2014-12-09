simulateNetwork<-function(n.segments, beta = NULL, 
                          lambda = 1, subs = F){  
  npts = 1
  network<-graph.tree(n.segments)
  coords<-layout.fruchterman.reingold.grid(network, circular = T)
  colnames(coords) <- c("easting", "northing")
  adjacency<-t(make_spam((get.adjacency(network))))
  wgts<-make_weights(adjacency)
  
  # eigenvalue decomp
  P<-spatial_penalty(adjacency, wgts, lambda = lambda, n.segments = n.segments)
  z<-eigen(P)
  zvals <- z$values
  zvecs <- z$vectors
  nzero <- which(zapsmall(zvals) == 0)
  nvals <- (length(zvals)-length(nzero))
  nr<-matrix(0, nrow = nvals, ncol = npts)
  for(k in 1:npts){
    nr[,k] <- rnorm(n=nvals, 
                mean = as.numeric(rep(0, nvals)), 
                sd = as.numeric(sqrt(1/zvals[-nzero]) )) 
  }
  z <- c(zvecs[,-nzero] %*% nr)
#   z<-z - mean(z)
  spatial <- z
  
  # add on linear covariate effects and intercept term
  if(length(beta) > 0){
    betaModel<-cbind(1, matrix(rnorm(n.segments*(length(beta) - 1)), 
           nrow = n.segments, ncol = (length(beta)-1)))
    betaContribution<-betaModel %*% beta
    z<-z + as.vector(betaContribution)
    betaNames<-"Intercept"
    if(!length(beta) == 1){
      for(i in 1:(length(beta)-1)){
        betaNames <- c(betaNames, paste("x", i, sep = ""))
      }
    }
    colnames(betaModel) <- betaNames
    betaModel<-Xout <-data.frame(betaModel)
    if(npts > 1){
      for(k in 1:(npts - 1)){
        Xout <- rbind(Xout, betaModel)
      }
    }
    betaModel <- Xout
  }
  
  # responses - linear part plus spatial network part plus random noise
  response <- as.vector(z + rnorm(n.segments*npts))
  response.locs = 1:(npts*n.segments)
  if(length(beta)>0){
    obsdata = data.frame(response, betaModel, coords) 
    } else  obsdata = data.frame(response, coords) 
  out<-list(obsdata = obsdata, response.locs = response.locs, 
            beta = beta, wgts  = wgts, adjacency = adjacency, spatial = spatial)
  
  # subset the data into a 'hold-out' prediction set and an observed set.  
  # Useful for testing the model prediction accuracy.
  if(!subs == F){
    subset.rows<-sample(x = 1:(n.segments*npts), size = round(n.segments*npts*subs, 0))
    Obs<-obsdata
    obsdata<-Obs[subset.rows,]
    preddata<-Obs[-(subset.rows),]
    all.response.locs<-response.locs
    response.locs<-all.response.locs[subset.rows]
    prediction.locs<-all.response.locs[-subset.rows]
#     npred<-length(prediction.locs)
#     prediction.matrix<-spam(list(i=1:npred, j=prediction.locs, rep(1, npred)), 
#                             nrow = npred, ncol = n.segments)
#     preddata<-cbind.spam(as.matrix(preddata), prediction.matrix)
    preddata<-data.frame(preddata, prediction.locs)
    adjacency<-adjacency
    out<-list(obsdata = obsdata, preddata = preddata, response.locs = response.locs, 
              prediction.locs = prediction.locs, spatial = spatial,
              wgts  = wgts, adjacency = adjacency, subs = subset.rows, coords = coords) 
  }
return(out)
}