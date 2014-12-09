
predict.smnet<-function(object, newdata = NULL, ...){
  # part of object that has spline smooth results
  smBit<-object[[2]]
   
  # if some prediction points are specified by the user
  if(!is.null(newdata)){
   out <- predictNewDataObject(object, newdata)
  } else {
    # if no points are spescified, use the SSN object
    out <- predictSSNobject(object)
  }
  return(out)
}