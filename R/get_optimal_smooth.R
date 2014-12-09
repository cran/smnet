
get_optimal_smooth<-function(P.list, X.spam, XTX.spam, X.list, 
                             response, control,
                             net = F, n.linear = n.linear,  
                             sm.names, lin.names, 
                             method, verbose, plot.pars = F, Pwee.list = Pwee.list){
  
  # SET UP COMPONENTS AND SIMPLE SUMMARIES REQUIRED FOR GETTING OPTIMAL PARAMETERS
  # ------------------------------------------------------------------------------
  # number of smoothing components (this is usually the same as or less than the number of smoothing pars)
  n.smooth.comps<-length(P.list)
  # vector where each element describes the number of smooth pars associated with each smooth component
  stretch_index<-vector("numeric", length = n.smooth.comps)
  for(i in 1:length(P.list)) stretch_index[i]<-ifelse(class(P.list[[i]]) == "list", length(P.list[[i]]), 1)
  # unnest the list of penalty matrices
  P.flat    <- make_flat(P.list)
  Pwee.flat <- make_flat(Pwee.list)
  # the total number of smoothing parameters (=sum(stretch_index))
  n.smooth.pars      <-length(P.flat)
  # length of the response vector
  n         <-length(response)
  # number of columns associated with each component of the data matrix
  X.dim     <- lapply(X.list, ncol)
  # number of components in the model - includes intercept and linear variables
  n.terms   <- length(X.list)
  # total columns of the data matrix
  np        <- sum(unlist(X.dim))
  # column numbers and cumulative column numbers in the data matrix
  inds      <- unlist(X.dim)
  cum.inds  <- cumsum(inds) 
  
  # IDENTIFIABILITY CONSTRAINTS USED IN ESTIMATION OF UNIVARIATE SMOOTHS
  # ------------------------------------------------------------------------------
  blockzero<-lapply(X.dim, make_sparse)
  idList <- XTX.list <- lapply(X.list, crossprodspam)
  makezero<-n.linear + 1
  # zero out the identifiability constraint associated with linear components, 
  #  and include those associated with standard spline smooths
  idList <- XTX.list
  ridgeM <- list()
  for(z in 1:makezero) idList[[z]] <- blockzero[[z]] 
  for(term in 1:(n.terms - net)) ridgeM[[term]]<-(10^-7)*get_block_penalty(P = idList[[term]], blockzero = blockzero, i=term)
  ridgeM <- Reduce("+", ridgeM)
  
  # IDENTIFIABILITY CONSTRAINTS USED FOR THE NETWORK COMPONENT, REGARDLESS OF "ridge".
  # THIS MATRIX IS FIXED THROUGHOUT, AND NO PARAMETER IS ASSOCIATED WITH IT.
  # ------------------------------------------------------------------------------
  nseg          <- inds[n.terms]
  net_par_start <- np - nseg + 1
  net_pars      <- net_par_start:np
  if(net)          diag.spam(ridgeM)[net_pars] <- 10^-7

  # CHOLESKY FACTORISATION, Xy.
  # ------------------------------------------------------------------------------
  P          <- P.flat
  Xy         <- t(X.spam) %*% response
  cholFactor <- chol.spam(XTX.spam + Reduce("+", P.flat) + ridgeM, pivot = "MMD", eps = 10^-6)
  
  # WORK OUT THE INITIAL VALUES OF THE SMOOTHING PARAMETERS FOR THE OPTIMISER
  # ------------------------------------------------------------------------------
  if(is.numeric(control$start)){
    start.vals <- control$start
  } else if(is.null(control$start)){
    start.vals <- rep(4, n.smooth.pars)
    if(net){
      start.vals[n.smooth.pars] <- -12
    }
  }

  if(method %in% c("AICC", "GCV", "AIC")){
      if(is.numeric(control$approx)){
        w <- matrix(sample(c(-1,1), n*control$approx, replace = T), nrow=n, ncol=control$approx)
        Xw<-t(X.spam) %*% w
        objective<-function(rhoArgs){
          get_crit_approx(rhoArgs, X = X.spam, 
                          XTX = XTX.spam, P = P, response = response,
                          Xy = Xy, Xw = Xw, cholFactor = cholFactor, n=n,
                          n.sm = n.smooth.pars, crit = method, identifyBit = ridgeM)
        }
      } else {
        objective<-function(rhoArgs){
          get_crit_exact(rhoArgs, X = X.spam, 
                         XTX = XTX.spam, P = P, response = response,
                         Xy = Xy, Xw = Xw, cholFactor = cholFactor, n=n,
                         n.sm = n.smooth.pars, crit = method, identifyBit = ridgeM)
        }
      }
      if(length(P) == 1){
        a1      <- proc.time()
        optimal <- optimize(f = objective, upper = 20, lower = -20, tol = control$tol)
        a2      <- proc.time()
        ntime   <- round(abs(a1-a2)[3], 3)
      }
      if(length(P) > 1){
        lower <- rep(-20, length(P.flat))
        upper <- rep(20, length(P.flat))
        if(!is.null(control$start)){
          if(control$start == "autostart"){
            if(length(P) == 2){
              n.pts<-5
              at.points.net<-seq(1, 14, length.out = n.pts)
              at.points.ridge<-seq(-14, 1, length.out = n.pts)
              z<-expand.grid(at.points.net, at.points.ridge)
              for(i in 1:nrow(z)){
                z[i,3]<-objective(rhoArgs = as.numeric(z[i,1:2]))
              }
              z<-as.matrix(z)
              mini<-which(z[,3] == min(z[,3], na.rm = TRUE))
              mini.best<-which(z[mini,1] == max(z[mini,1]))
              mini.grid.pars<-as.numeric(z[mini[mini.best],1:2])
              if(length(mini.grid.pars)>2) mini.grid.pars <- mini.grid.pars[c(1,3)]
              start.vals<-mini.grid.pars
            }
          }
        }
        if(net){
          lower[n.smooth.pars - 1] <- 0
          upper[n.smooth.pars - 1] <- 30
          lower[n.smooth.pars] <- -20
          upper[n.smooth.pars] <- 2
        }
        a1<-proc.time()
        optimal<-nmkb(par = start.vals, fn = objective, 
                      lower = lower, upper = upper, control=list(tol = control$tol))
        a2<-proc.time()
        ntime<-round(abs(a1-a2)[3], 3)
        nits<-ifelse(length(P) > 1, optimal$feval, 1)
      } 
  }
  
  # RETRIEVE THE OPTIMAL SMOOTHING PARAMETER VALUES FROM THE ABOVE
  # ------------------------------------------------------------------------------
  if(length(P.flat) == 1) rhoOptim <- optimal$minimum
  if(length(P.flat) > 1)  rhoOptim <- optimal$par
  
  # FINALLY FIT THE MODEL USING THE ESTIMATED SMOOTHING PARAMETERS
  # ------------------------------------------------------------------------------
  for(j in 1:n.smooth.pars) P[[j]]<-P[[j]]*exp(rhoOptim[j])
  Psum     <- Reduce("+", P)
  info     <- XTX.spam + Psum + ridgeM
  U        <- update.spam.chol.NgPeyton(cholFactor, info)
  beta_hat <- backsolve.spam(U, forwardsolve.spam(U, Xy))  
  left1    <- forwardsolve.spam(U, t(X.spam))
  trH      <- rowSums(left1*left1)
  ED1      <- sum(trH)
  left2    <- backsolve.spam(U, left1)
  left2    <- X.spam %*% left2
  ED2      <- sum(left2*left2)
  dof.list <- vector("list", length = n.terms)
  dof.sm   <- vector("list", length = n.smooth.comps)
  dof.list[[1]]<-trH[1]
  for(i in 2:n.terms)           dof.list[[i]] <- trH[(cum.inds[i-1]+1):(cum.inds[i])]
  for(i in 1:n.smooth.comps)    dof.sm[[i]]   <- dof.list[[which(!inds == 1)[i]]]
  fit      <- X.spam %*% beta_hat
  
  # VARIANCE AND DIAGNOSTICS FOR OUTPUT
  # ------------------------------------------------------------------------------
  sigma.sq    <- sum((response - fit)^2)/(n - 2*ED1 + ED2) # model variance
  pdof        <- lapply(dof.sm, sum)
  dof.smooth  <- rep(pdof, stretch_index)
  allrho      <- exp(rhoOptim)
  diagnostics <- list(ED = 2*ED1 - ED2, sigma.sq = sigma.sq, fit = fit)
  
  # IF VERBOSE IS SET TO TRUE, PRINT THE CONVERGENCE STATUS FROM OPTIMISER
  # ------------------------------------------------------------------------------
  if(length(P.flat) > 1){
    if(control$verbose){
      if(optimal$convergence == 0){
        nits<-ifelse(length(P) == 1, optimal$counts, optimal$feval)
        print(paste("Convergence reached in", nits, "iterations after ", ntime, " seconds"))
      } else { print ("Warning, convergence was not reached")}   
    }
  } else if(length(P.flat) == 1) print(paste("Convergence reached after ", ntime, " seconds"))
  
  # PRINT THE SMOOTHING PARAMETERS, AND ASSOCIATED EFFECTIVE DIMENSIONS
  # ------------------------------------------------------------------------------
  get.inner.length<-function(L) ifelse(!class(L) == "list", 1, length(L))
  inner.length<-lapply(P.list, get.inner.length)
  tablecol<-max(unlist(inner.length))
  tablerow<-length(sm.names)+net
  retPrint<-ret<-matrix(as.factor("---"), nrow = tablerow, ncol = tablecol)
  increment<-1
  for(i in 1:length(P.list)){
    for(j in 1:inner.length[[i]]){
      if(round(allrho[increment], 3) > 0) nice_number <- round(allrho[increment], 3)
      if(10^-8 >  allrho[increment])                             nice_number <- "<10^-8"
      if(10^-8 <= allrho[increment] & allrho[increment] < 10^-7) nice_number <- "~10^-7"
      if(10^-7 <= allrho[increment] & allrho[increment] < 10^-6) nice_number <- "~10^-6"
      if(10^-6 <= allrho[increment] & allrho[increment] < 10^-5) nice_number <- "~10^-5"
      if(10^-5 <= allrho[increment] & allrho[increment] < 10^-4) nice_number <- "~10^-4"
      if(10^8  <  allrho[increment])                             nice_number <- ">10^8"
      retPrint[i,j]  <- nice_number
      ret[i,j]       <- allrho[increment]
      increment      <- increment+1
    }
  }
  smpar.names<-vector(length = tablecol)
  for(i in 1:tablecol) smpar.names[i]<-paste("Lambda_", i, sep = "")
  retPrint            <- cbind(retPrint, round(unlist(pdof), 2))
  ret                 <- cbind(ret, round(unlist(pdof), 2))
  rownames(retPrint)  <- c(if(net){c(sm.names, "Network")}else{sm.names})
  rownames(ret)       <- c(if(net){c(sm.names, "Network")}else{sm.names})
  colnames(retPrint)  <- c(smpar.names, "pED")
  colnames(ret)       <- c(smpar.names, "pED")
  if(control$verbose) print(as.data.frame(retPrint)) 
  
  # RETURN MODEL OUTPUT
  # ------------------------------------------------------------------------------
  list(pars = allrho, beta_hat = beta_hat, ED = 2*ED1 - ED2, fit = fit, 
       out = ret, sigma.sq = diagnostics$sigma.sq, U = U)
}


