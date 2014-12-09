
plot_node<-function(fitted_segment, adjacency, node_coords, 
                    xlim=range(node_coords[,1]), ylim=range(node_coords[,2]), 
                    col = 3, cex = NULL, key = TRUE, log = NULL)
{
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
#   hist(fitted_segment
  nnodes <- length(fitted_segment)
  nclasses <- max(c(10, round(nnodes/10, 0)))
  lower <- min(fitted_segment, na.rm = T) - .0001
  upper <- max(fitted_segment, na.rm = T) + .0001
  brks <- seq(lower, upper, length.out = nclasses)
#   brks <- quantile(fitted_segment, probs = (1:(nclasses-1))/nclasses, na.rm = T)

  colorPalette<-heat.colors(nclasses-1)
  col.nums<-cut(as.vector(fitted_segment), breaks = brks, labels = FALSE)
  coord = node_coords
  col = colorPalette[col.nums]
  if(is.null(cex)) cex = 020/sqrt(nnodes)
  layout(matrix(1:2, nrow = 1), widths = c(4,1))
  par(mar = c(0,0,0,1))
  plot(coord[,2]~coord[,1], type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty="n", xlim = xlim, ylim=ylim)
  for(i in 1:ncol(adjacency)){
    nums<-which(!adjacency[,i] == 0)
    if(length(nums)>0){
      for(j in 1:length(nums)){
        segments(x0=coord[i,1], y0=coord[i,2],
                 x1=coord[nums[j],1], y1=coord[nums[j],2], lwd = 1)
      }
    }
  }
  points(coord, pch = 21, bg =  "black", cex = as.vector(cex)+0.1, col = 1)
  points(coord, pch = 21, bg =  as.vector(col), cex = as.vector(cex), col = 1)
  par(mar = c(1, 1.5, 1, 3))
  if(key) colour.key(brks = brks, cols = colorPalette)
par(def.par)
}