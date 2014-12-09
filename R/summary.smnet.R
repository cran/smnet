summary.smnet<-function(object,...){
  smPart<-object[[2]]
  adjacency<-smPart$adjacency
  n.smooth<-smPart$n.smooth
  out<-smPart$out
  n.linear<-smPart$n.linear
  X.list<-smPart$X.list
  beta_hat<-smPart$beta_hat
  U<-smPart$U
  sigma.sq<-smPart$sigma.sq
  XTX.spam<-smPart$XTX.spam
  X.spam<-smPart$X.spam
  lin.names<-smPart$lin.names
  ED<-smPart$ED
  fit<-smPart$fit
  n <- length(smPart$response)
  response <- smPart$response
  
  # Summarise the smooth part of the model
  smoothExists<-n.smooth + !is.null(adjacency)
  outList<-vector("list")
  if(smoothExists>0) outList$smoothcomponents<-as.data.frame(out)
  
  # Summarise the linear part of the model
  linearExists<-n.linear + 1
  if(linearExists>0){
    getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
    cov.dims  <- lapply(X.list, getcol)    
    inds      <- unlist(cov.dims)
    cum.inds  <- cumsum(inds)
    n.cov     <- length(X.list)
    unit.locations<-which(inds == 1)
    
    # get standard errors 
    left1<-forwardsolve.spam(U, t(X.spam))
    left2<-backsolve.spam(U, left1)
    separs<-sqrt(sigma.sq*rowSums(left2^2)[unit.locations])
    pars<-beta_hat[unit.locations]
    tval<-round(pars/separs, 6)
    pars<-round(pars, 6)
    separs<-round(separs, 6)
    ppar <- (1 - pt(abs(pars/separs), df = n - ED)) * 2
    R2 <- 1 - (sum((response - fit)^2)/sum((response - mean(response))^2))
    print(R2)
    tablecol<-4
    tablerow<-linearExists
    ret<-matrix(0, nrow = tablerow, ncol = tablecol)
    increment<-1
    rownames(ret)<-c("Intercept", lin.names) 
    colnames(ret)<-c("Estimate", "Std.Err", "t-Value", "p-value")
    ret[,1]<-pars
    ret[,2]<-separs
    ret[,3]<-tval
    ret[,4]<-ppar
    outList$fixed.effects<-as.data.frame(ret)
    return(outList)
  }
}




  