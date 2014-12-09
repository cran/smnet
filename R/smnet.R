
smnet<-function(formula, data.object, control = NULL, method = "AICC")
{  
  
  adjacency <- NULL
  # extract data if SSN data is provided
  if(class(data.object) == "SpatialStreamNetwork"){
    data            <- getSSNdata.frame(data.object, Name = "Obs")
  } 
  # if data is provided as a dataframe
  if(class(data.object) == "data.frame"){
    data <- data.object
  }
  # error-check the formula, extracting all data to construct the smooths
  formulaout<-try(formula_stuff <- get_formula_stuff(formula = formula, data = data), silent=T)
  if(class(formulaout) == "try-error"){
    variable.names<-names(data)
    cat("One or more variables does not exist, check they are all one of the following:")
    print(variable.names)
    stop()
  } 
  
  # construct the vector of weights, apply an ordering if appropriate
  cls.sm<-lapply(formulaout$gp$smooth.spec, class)
  net <- "network.spec" %in% cls.sm
  if(net){
    networkSpecNo<-which(cls.sm== "network.spec")
    networkSpec<-formulaout$gp$smooth.spec[[networkSpecNo]]
    adjacency<-networkSpec$adjacency
    if(!is.null(adjacency)){
      if(class(data.object) == "SpatialStreamNetwork"){
        netID<-networkSpec$netID
        # get the rids out of the SSN object associated with netID
        data            <- data[as.numeric(as.character(data$netID)) == netID,]
        rid_data        <- data.object@data
        rid_data        <- rid_data[rid_data$netID == netID,]
        if(nrow(rid_data) == 0) stop("No data associated with that netID")
        ord             <- order(rid_data$rid)
        response.locs   <- as.numeric(as.character(data$rid))
        # if the weights are unprocessed stream order like shreve...
        weight <- networkSpec$weight
        if(weight %in% c("shreve", "addfunccol")){
          weight       <- try(rid_data[weight], silent = T)
          if(class(weight) == "try-error") stop("No data found for specified weight")
          weight       <- get_shreve_weights(adjacency = adjacency, shreve.order = as.matrix(weight))
        } else if(!(weight %in% c("shreve", "addfunccol"))){
          # if the weights are already processed like areaPI...
          weight         <- try(as.vector(as.matrix(rid_data[weight])), silent = T)
          if(class(weight) == "try-error") stop("No data found for specified weight")
          weight         <- weight[ord]
        } 
      }
      if(class(data.object) == "data.frame"){
        if(is.character(networkSpec$locs)){
          response.locs <- try(data[networkSpec$locs], silent = T)
          if(response.locs == "try-error"){
            stop("No data found for locs provided")
          }
        } else if(is.numeric(networkSpec$locs)){
          response.locs <- networkSpec$locs
        } else stop("locs should be numeric vector or a character string")
        weight <- networkSpec$weight
        if(!is.numeric(weight)) stop("If no data is not an SSn obejct, weight must be a vector")
        ord <- rid_data <- netID <- NULL
      }
      # check that weights vector and adjacency are the same dimension
      if(!nrow(adjacency) == length(weight)) stop("weights and adjacency have different dimensions")
    } else ord <- rid_data <- netID <- NULL
  } else ord <- rid_data <- netID <- weight <- NULL


  default.control = list(trace = 0, maxit = 500, start = NULL, approx = NULL, verbose = TRUE, tol = 10^-5)
  if(!is.null(control)) for(i in 1:length(control)) default.control[names(control)[i]]<-list(control[[i]])
  # create model objects
  model_objects <- get_model_objects(formula = formula, data = data, 
                                     adjacency = adjacency, 
                                     response.locs = response.locs, weight = weight,  
                                     rid_data = rid_data, netID = netID, ord = ord, 
                                     control = default.control, formulaout = formulaout)
  
  #  choose optimal smooth parameters using box constrained Nelder-Mead search
  opt <- with(model_objects, {
    if((n.smooth == 0) && (!net)){
      get_lm_fit(X.spam = X.spam, X.list = X.list, XTX.spam = XTX.spam, 
                 response = response)         
    } else {
      get_optimal_smooth(P.list = P.list, X.spam = X.spam, X.list=X.list, XTX.spam=XTX.spam, 
                         response=response, control = control, net=net, n.linear = n.linear,
                         lin.names = lin.names, sm.names = sm.names, 
                         method = method, Pwee.list = Pwee.list, plot.pars = plot.pars)   
    }
#       P.list<-model_objects$P.list
#     Pwee.list<-model_objects$Pwee.list
#       X.spam<-model_objects$X.spam
#       X.list<-model_objects$X.list
#       XTX.spam<-model_objects$XTX.spam
#       response<-model_objects$response
#       control<-model_objects$control
#         n.linear<-model_objects$n.linear
#       sm.names<-model_objects$sm.names
#       smp<-model_objects$smp
#       crit<-model_objects$crit
#       verbose<-model_objects$verbose

  }
  )
  
  # create output list
  outputList<-vector("list", length = 2)
  outputList[[1]]<-data.object
  output<-c(opt, model_objects)
 

  outputList[[2]]<-output
  class(outputList)<-"smnet"
  outputList
}

