
plot_continuous<-function(x, fitted_segment, ord){  	
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  layout(matrix(1:2, nrow = 1), widths = c(4,1))
  par(mar = c(5,5,3,0))
  plot(x@bbox[1,],x@bbox[2,], type = "n")
  data <- x@data
  # create quantile breaks
  reorder<-order(ord)
  nclasses<-20
  brks <- quantile(fitted_segment,
                   probs = (1:(nclasses-1))/nclasses, na.rm = T)
  lower.breaks <- c(min(fitted_segment, na.rm = T) - .0001, brks)
  upper.breaks <- c(brks, max(fitted_segment, na.rm = T))
  colorPalette <- rainbow(nclasses, start = .66, end = .99)
  for(k in 1:nclasses) {
    for(i in 1:length(x@lines))
      for(j in 1:length(x@lines[[reorder[i]]]))
        if(fitted_segment[i] > lower.breaks[k] & fitted_segment[i] <= upper.breaks[k])
          lines(x@lines[[reorder[i]]]@Lines[[j]]@coords,
                col = colorPalette[k])
  }
  dec.dig <- 2
  left <- as.character(as.numeric(as.integer(lower.breaks*10^dec.dig))/
                         10^dec.dig)
  rght <- as.character(as.numeric(as.integer(upper.breaks*10^dec.dig))/
                         10^dec.dig)
  leglabs <- paste(left,"to",rght)
  par(mar = c(0,0,0,0))
  plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab ="", bty = "n")
  legend(x = -1, y = 1.1, legend = leglabs, bty = "n",
         lty = rep(1, times = nclasses),
         col = colorPalette, cex = .8)
  par(def.par)
}

