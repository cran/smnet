plot.smnet<- function(x, 
                      type = "covariates", 
                      se = FALSE, 
                      res = FALSE, 
                      coords = NULL,
                      key = TRUE, ...)
{
  if(!class(x) == "smnet") stop("x must be an object of class 'smnet'")
  if((type %in% c("nodes", "segments")) && (res == T)) warning("Ignoring 'res' argument, since spatial plot requested")
  if(type == "nodes"){
    if(!x[[2]]$net) stop("No spatial network component to plot in x")
      # break up beta_hat
      with(x[[2]], {
        getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
        cov.dims  <- lapply(X.list, getcol)    
        inds      <- unlist(cov.dims)
        cum.inds  <- cumsum(inds)
        n.cov     <- length(X.list)
        ests      <- vector("list", length = n.cov)
        ests[[1]] <- beta_hat[1]
        for(i in 1:(n.cov-1)){ests[[i+1]]<-beta_hat[(cum.inds[i]+1):(cum.inds[i+1])]} 
        spatial.comp<-ests[[n.cov]]# + ests[[1]]
        if(is.null(coords)){
          # get midpoint locations for each rid, then put in order
          data      <- getSSNdata.frame(x[[1]], Name = "Obs")
          get_rid_midpt<-function(L) apply(L@Lines[[1]]@coords, 2, FUN = "mean") 
          midpt_rids<-lapply(x[[1]]@lines, FUN = get_rid_midpt)
          coords<-Reduce("rbind", midpt_rids)[ord,]
        } 
        plot_node(fitted_segment=spatial.comp, adjacency = adjacency, 
                  node_coords=coords, key = key,...)
      })
  }
  if(type == "segments"){
    if(!x[[2]]$net) stop("No spatial network component to plot in x")
    # break up beta_hat
    with(x[[2]], {
      getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
      cov.dims  <- lapply(X.list, getcol)    
      inds      <- unlist(cov.dims)
      cum.inds  <- cumsum(inds)
      n.cov     <- length(X.list)
      ests      <- vector("list", length = n.cov)
      ests[[1]] <- beta_hat[1]
      for(i in 1:(n.cov-1)) ests[[i+1]] <- beta_hat[(cum.inds[i]+1):(cum.inds[i+1])]
      spatial.comp <- ests[[n.cov]]# + ests[[1]]
      plot_continuous(x[[1]], fitted_segment=spatial.comp, ord=ord) 
    })
  }
  if(type == "covariates"){
    plot_cov(x=x[[1]], spline_part = x[[2]], res = res, ...)
  }
}
 


