get_crit_approx<-function(rhoArgs, X, XTX, P, response, cholFactor, n, 
                          np = nrow(XTX), Xw, Xy, n.sm, identifyBit, crit){    
  for(j in 1:n.sm) P[[j]]<-P[[j]]*exp(rhoArgs[j])
  Psum<-Reduce("+", P)
  info<-XTX+Psum+identifyBit
  U<-update.spam.chol.NgPeyton(cholFactor, info) 
  beta_hat<-backsolve.spam(U, forwardsolve.spam(U, Xy)) 
  left1<-forwardsolve.spam(U, Xw)
  left2<-backsolve.spam(U, left1)
  left2<-X %*% left2
  ED<-sum(left2*left2)
  trH<-as.vector(apply(as.matrix(left1*left1), 1, "mean"))
  ED<-2*sum(trH) - ED
  if(ED > n){n <- ED + 1}
  ED<-sum(trH)
  fit<-X %*% beta_hat  
  if(crit == "AICC") {
    out<-log(mean((response - fit)^2)) + (2*(ED+1)/(n-ED))
  }
  else if(crit == "GCV"){
    out <- mean((response - fit)^2)/(1-ED/n)^2
  } else if(crit == "AIC"){
    out<-log(mean((response - fit)^2)) + 2*ED/n
  }
  out
}

# 


get_ED_approx<-function(rhoArgs, X, XTX, P, response, cholFactor, n, 
                 np = nrow(XTX), Xw, Xy, n.sm, identifyBit, crit){    
  for(j in 1:n.sm) P[[j]]<-P[[j]]*exp(rhoArgs[j])
  Psum<-Reduce("+", P)
  info<-XTX+Psum+identifyBit
  U<-update.spam.chol.NgPeyton(cholFactor, info) 
  beta_hat<-backsolve.spam(U, forwardsolve.spam(U, Xy)) 
  left1<-forwardsolve.spam(U, Xw)
    left2<-backsolve.spam(U, left1)
    left2<-X %*% left2
    ED<-sum(left2*left2)
  trH<-as.vector(apply(left1*left1, 1, "mean"))
    ED<-2*sum(trH) - ED
    if(ED > n){n <- ED + 1}
  sum(trH)
}