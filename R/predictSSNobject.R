
predictSSNobject<-function(smNetFit){
  #   seperate the SSN object and the spline model fit
  SSNPart<-smNetFit[[1]]
  smPart<-smNetFit[[2]]
  
  #   put all of the components of new design matrix in a list
  newX<-vector("list")
  
  # extract the prediction locations and the prediction data matrix
  if(class(SSNPart) == "SpatialStreamNetwork"){
    rids<-SSNPart@data$rid#[SSNPart@data$netID == netID]
    rid_ord<-rids[order(rids)]
    dfPred<-getSSNdata.frame(SSNPart, "preds")
    ridPred<-as.numeric(as.character(dfPred[,"rid"]))
  }
  
  # get the linear components if it exists
  if(smPart$n.linear > 0){
    linVarbs<-dfPred[,colnames(smPart$variables)[1:smPart$n.linear]]
    newX<-c(newX,  list(linear = cbind(1, linVarbs)))
  } else {
    newX<-c(newX, list(linear = rep(1, length(ridPred))))
  }
  
  #   get the smooth additive components if there are any  
  if(smNetFit[[2]]$n.smooth > 0){
    smTerms<-smNetFit[[2]]$sm.terms.names
    smDesign<-vector("list")
    dataOriginal<-data.frame(smNetFit[[2]]$variables)
    for(i in 1:length(smTerms)){
      oldVariable<-dataOriginal[smTerms[[i]]]
      if(length(smTerms[[i]]) == 1){
        xlxr1<-range(oldVariable)
        newX<-c(newX, b_spline_basis(x=unlist(dfPred[smTerms[[i]]]), 
                                     nseg = (smNetFit[[2]]$sm.basis[i] - 3), deg=3))
      }
      if(length(smTerms[[i]]) == 2){
        xlxr1<-range(oldVariable[,1])
        xlxr2<-range(oldVariable[,2])
        a1<-b_spline_basis(x=unlist(dfPred[smTerms[[i]][1]]), 
                           xl=xlxr1[1], xr=xlxr1[2],nseg = (smNetFit[[2]]$sm.basis[i]-3), deg=3)
        a2<-b_spline_basis(x=unlist(dfPred[smTerms[[i]][2]]), 
                           xl=xlxr2[1], xr=xlxr2[2],nseg = (smNetFit[[2]]$sm.basis[i]-3), deg=3)
        newX<-c(newX, make_spam(not_sparse_box_product(a1, a2)))
      }
    }
  }
  
  if(smPart$net){
    # construct network component
    newX<-c(newX, spam(x = list(i=1:length(ridPred), j=ridPred, 
                                val = rep(1, length(ridPred))), nrow = length(ridPred), 
                       ncol = ncol(smPart$X.list[[length(smPart$X.list)]])))
  }
  
  #   put together the linear, smooth and network components
  newX<-lapply(newX, as.matrix)
  Xstar<-Reduce("cbind.spam", newX)
  predictions<-Xstar %*% smPart$beta_hat
  left1<-forwardsolve.spam(smPart$U, t(smPart$X.spam))
  left2<-backsolve.spam(smPart$U, left1)
  vec<-  Xstar %*% left2
  predictions.se<-sqrt((1 + rowSums(vec*vec))*smPart$diagnostics$sigma.sq)
  list(predictions = predictions, predictions.se = predictions.se)   
  }
