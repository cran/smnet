plot_cov<-function(x, spline_part, res = res, ...)
{
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  with(spline_part, 
{
  X.list<-spline_part$X.list
  beta_hat<-spline_part$beta_hat
  variables<-spline_part$variables
  sigma.sq    <- spline_part$sigma.sq
  n.linear<-spline_part$n.linear
  n.smooth<-spline_part$n.smooth
  varbs.len<-spline_part$varbs.len
  varbs.loc<-spline_part$varbs.loc
  sm.cyclic<-spline_part$sm.cyclic
  sm.basis<-spline_part$sm.basis
  U<-spline_part$U
  X.spam<-spline_part$X.spam
  response <- spline_part$response
  fit      <- spline_part$fit
  getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
  cov.dims  <- lapply(  X.list, getcol) 
  np        <- sum(unlist(cov.dims))
  cs.np     <- cumsum(unlist(cov.dims))
  n.cov     <- length(X.list)
  n.segments<- cov.dims[[n.cov]]
  
  # break up beta_hat
  inds      <- unlist(cov.dims)
  cum.inds  <- cumsum(inds)
  ests      <- vector("list", length = n.cov)
  ests[[1]]<-beta_hat[1]
  for(i in 1:(n.cov-1)){ests[[i+1]]<-beta_hat[(cum.inds[i]+1):(cum.inds[i+1])]}
  cov.names   <- colnames(variables)
  n.varbs     <- ncol(variables)
  fit.mean    <- vector("list", length = n.varbs)
  partial_fit <- vector("list", length = n.varbs)
  new.varb    <- vector("list", length = n.varbs)
  se.list     <- se.ords<-vector("list", length = n.varbs)
  
  # create standard errors on regular grid; linear terms first
  if(n.linear>0){
    for(i in 1:n.linear){
      varb<-variables[,i]
      block.inds <- c((cs.np[i]+1):cs.np[(i+1)])
      new.X<-create_component_design(varb=varb,X.list=X.list,j=(i+1),sm=F)
      new.X[,1]<-1
      #new.X      <- new.X[,block.inds]
      left1<-forwardsolve.spam(U, t(new.X), transpose = T)
      left2<-backsolve.spam(U, left1, transpose = T)
      vec<- X.spam %*% left2
      se.list[[i]]<-sqrt(colSums(vec*vec)*sigma.sq)
    }
  }
  if(n.smooth>0){
    new.basis<-vector("list", length = n.smooth)
    new.varb.smooth<-vector("list", length = n.smooth)
    for(i in 1:n.smooth){
      if(varbs.len[i] == 1){
        varb       <- variables[,varbs.loc[[i]]]
        block.inds <- c( (cs.np[i+n.linear]+1):cs.np[(i+n.linear+1)])
        new.varb.smooth[[i]]<-seq(min(varb),  max(varb), length.out = 100)
        new.basis[[i]]<-make_spam(bbase(new.varb.smooth[[i]], nseg = (sm.basis[i]-3), deg = 3))
        new.X<-create_component_design(varb=new.varb.smooth[[i]],X.list=X.list,j=(i+1+n.linear), sm=T)
        left1<-forwardsolve.spam(U, t(new.X), transpose = T)
        left2<-backsolve.spam(U, left1, transpose = T)
        vec<- X.spam %*% left2
        se.list[[i+n.linear]]<-sqrt(colSums(vec*vec)*sigma.sq)   
      }
      if(varbs.len[i] == 2){
        locs       <- varbs.loc[[i]]
        varb       <- variables[,locs]
        new.varb.1 <- seq(min(varb[,1]),  max(varb[,1]), length.out = 50)
        new.varb.2 <- seq(min(varb[,2]),  max(varb[,2]), length.out = 50)
        new.bas.1  <- make_spam(bbase(new.varb.1, nseg = (sm.basis[i]-3), deg = 3))
        new.bas.2  <- make_spam(bbase(new.varb.2, nseg = (sm.basis[i]-3), deg = 3))
        new.bas.1  <- new.bas.1 %x% matrix(1, nrow = 50, ncol = 1)
        new.bas.2  <- matrix(1, nrow = 50, ncol = 1) %x% new.bas.2
        new.X      <- new.basis[[i]] <- make_spam(not_sparse_box_product(new.bas.1, new.bas.2))
        X.dim <- lapply(X.list, ncol)
        empty.X.list <- lapply(X.dim, make_sparse, nrow = 2500)
        empty.X.list[[i+n.linear+1]] <- new.X
        new.X<-Reduce("cbind", empty.X.list)
        left1<-forwardsolve.spam(U, t(new.X), transpose = T)
        left2<-backsolve.spam(U, left1, transpose = T)
        vec<- X.spam %*% left2
        se.list[[i+n.linear]]<-sqrt(colSums(vec*vec)*sigma.sq)   
#         block.inds <- c((cs.np[i+n.linear]+1):cs.np[(i+n.linear+1)])
#         se.list[[i+n.linear]]<-sqrt(diag.spam(new.X %*% var.mat[block.inds,block.inds] %*% t(new.X)) * sigma.sq)
      }
    }
  }
  # create partial fits on regular grid; linear terms first
  if(n.linear>0){
    for(i in 1:n.linear){
      varb<-variables[,i]
      rnge<-range(varb)
      new.varb[[i]]<-as.matrix(seq(rnge[1], rnge[2], length.out = 100))
      fit.mean[[i]]<-new.varb[[i]]%*%ests[[i+1]] + ests[[1]]
      partial_fit[[i]]<-varb*ests[[i+1]] + ests[[1]]
  
    }
  }
  if(n.smooth){
    for(i in 1:n.smooth){
      fit.mean[[i+n.linear]]<-new.basis[[i]]%*%ests[[i+1+n.linear]]
      partial_fit[[i+n.linear]]<-X.list[[i+1+n.linear]]%*%ests[[i+1+n.linear]]
    }
  } 
  par(...)
  # get max and min for plots, plotting effects on same scale
  mn<-0
  mx<-0
  if(n.linear>0){
    for(i in 1:n.linear){
      upper<-2*se.list[[i]]+fit.mean[[i]]
      lower<--2*se.list[[i]]+fit.mean[[i]]
      mn<-ifelse(min(lower)<mn, min(lower), mn)
      mx<-ifelse(max(upper)>mx, max(upper), mx)
    }
  }
  if(n.smooth>0){
    for(i in 1:n.smooth){
      upper<- 2*se.list[[i+n.linear]]+fit.mean[[i+n.linear]]
      lower<- -2*se.list[[i+n.linear]]+fit.mean[[i+n.linear]]
      mn<-ifelse(min(lower)<mn, min(lower), mn)
      mx<-ifelse(max(upper)>mx, max(upper), mx)
    }
  }
  
  # Plot linear first
  if(n.linear>0){
    for(i in 1:n.linear){
      partial_resids<-(response-fit) + partial_fit[[i]]
      type<-ifelse(!res, "n", "p")
      mn<-min(min(-2*se.list[[i]]+fit.mean[[i]]), min(partial_resids))
      mx<-max(max(2*se.list[[i]]+fit.mean[[i]]), max(partial_resids))
      plot.default(variables[,i], partial_resids, ylim=c(mn, mx), type = type, 
                   xlab = cov.names[i], ylab = "")  
      lines(new.varb[[i]], 2*se.list[[i]]+fit.mean[[i]], lty=3, lwd = 1.5, col=3)
      lines(new.varb[[i]],-2*se.list[[i]]+fit.mean[[i]], lty=3, lwd = 1.5, col =3)
      lines(new.varb[[i]], fit.mean[[i]], col = "red", lwd = 2)
      if(!res) rug(variables[,i])
    }
  }
  if(n.smooth>0){
    for(i in 1:n.smooth){
      if(varbs.len[i] == 1){
        partial_resids<-(response-fit) + partial_fit[[i+n.linear]]
        upperCI<- 2*se.list[[i+n.linear]]+fit.mean[[i+n.linear]]
        lowerCI<- -2*se.list[[i+n.linear]]+fit.mean[[i+n.linear]]
        ylim.v<-range(c(upperCI, lowerCI))
        type<-ifelse(!res, "n", "p")
        mn<-min(mn, min(partial_resids))
        mx<-max(mn, max(partial_resids))
        plot.default(variables[,i+n.linear], partial_resids, type = type, ylim=c(mn, mx), 
                     xlab = cov.names[i+n.linear], ylab = "")  
        lines(new.varb.smooth[[i]], upperCI, lty=2, lwd = 1, col=1)
        lines(new.varb.smooth[[i]], lowerCI, lty=2, lwd = 1, col =1)
        lines(new.varb.smooth[[i]], fit.mean[[i+n.linear]], col = "red", lwd = 2)
        if(!res) rug(variables[,i+n.linear])
      }
      if(varbs.len[i] == 2){
        # plotting code for 2d smooth interaction terms, mean and standard errors
        partial_resids<-(response-fit) + partial_fit[[i+n.linear]]
        ylim.v<-range(c(partial_resids, fit.mean[[i+n.linear]]-2*se.list[[i+n.linear]],
                        fit.mean[[i+n.linear]]-2*se.list[[i+n.linear]]))
        varbs<-variables[,varbs.loc[[i]]]
        index.1<-seq(min(varbs[,1]), max(varbs[,1]), length.out = 50)
        index.2<-seq(min(varbs[,2]), max(varbs[,2]), length.out = 50)
        i.n<-cov.names[varbs.loc[[i]]]
        z<-matrix(fit.mean[[i+n.linear]], nrow = 50, ncol = 50, byrow = T)
        z.se<-matrix(se.list[[i+n.linear]], nrow = 50, ncol = 50, byrow = T)
        cols<-heat.colors(15)
        brksm<-seq(min(z), max(z), length.out = (length(cols)+1))
        brkse<-seq(min(z.se), max(z.se), length.out = (length(cols)+1))
        layout(matrix(1:2, nrow = 1), widths = c(4,1))
        image(x = index.1, y = index.2, z=z, main = "Mean",xlab = i.n[1],ylab=i.n[2]) #col=cols,breaks=brksm)
        par(mar = c(2,4,2,3))
        colour.key(brks = brksm, cols = cols)
        layout(matrix(1:2, nrow = 1), widths = c(4,1))
        image(x = index.1, y = index.2, z=z.se, main = "Standard errors",xlab = i.n[1],ylab=i.n[2])#, col=cols,breaks=brkse)
        par(mar = c(2,4,2,3))
        colour.key(brks = brkse, cols = cols)
      }
    }
  }
})
par(def.par)
}