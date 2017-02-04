SpectrumSimilarity <- function(spec.top, spec.bottom, t = 0.25, b = 10, top.label = NULL, 
                               bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, print.alignment = FALSE,
                               print.graphic = TRUE, output.list = FALSE) {
  
  
  ## format spectra and normalize intensitites
  
  top_tmp <- data.frame(mz = spec.top[,1], intensity = spec.top[,2])
  top_tmp$normalized <- round((top_tmp$intensity / max(top_tmp$intensity)) * 100)
  top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= xlim[2])   
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)   # data frame for plotting spectrum
  top <- subset(top_plot, top_plot$intensity >= b)   # data frame for similarity score calculation
  
  bottom_tmp <- data.frame(mz = spec.bottom[,1], intensity = spec.bottom[,2])
  bottom_tmp$normalized <- round((bottom_tmp$intensity / max(bottom_tmp$intensity)) * 100)
  bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & bottom_tmp$mz <= xlim[2])   
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)   # data frame for plotting spectrum
  bottom <- subset(bottom_plot, bottom_plot$intensity >= b)   # data frame for similarity score calculation
  
  
  ## align the m/z axis of the two spectra, the bottom spectrum is used as the reference
  
  for(i in 1:nrow(bottom))
    top[,1][bottom[,1][i] >= top[,1] - t & bottom[,1][i] <= top[,1] + t] <- bottom[,1][i]
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  if(length(unique(alignment[,1])) != length(alignment[,1])) warning("the m/z tolerance is set too high")
  alignment[,c(2,3)][is.na(alignment[,c(2,3)])] <- 0   # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  
  if(print.alignment == TRUE) {
    print(alignment)
  }
  
  
  ## similarity score calculation
  
  if(x.threshold < 0) stop("x.threshold argument must be zero or a positive number")
  
  alignment <- alignment[alignment[,1] >= x.threshold, ]
  
  u <- alignment[,2]; v <- alignment[,3]
  similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
  
  
  ## generate plot
  
  if(print.graphic == TRUE) {
    plot.new()
    plot.window(xlim = xlim, ylim = c(-125, 125))
    ticks <- c(-100, -50, 0, 50, 100)
    for(i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]), col = "blue")
    for(i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]), col = "red")
    axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
    axis(1, pos = -125)
    lines(xlim, c(0,0))
    rect(xlim[1], -125, xlim[2], 125)
    mtext("m/z", side = 1, line = 2)
    mtext("intensity (%)", side = 2, line = 2)
    plot.window(xlim = c(0, 20), ylim = c(-10, 10))
    text(10, 9, top.label)
    text(10, -9, bottom.label)
  }
  
  if(output.list == TRUE) {
    
    # Grid graphics head to tail plot
    
    headTailPlot <- function() {
      
      pushViewport(plotViewport(c(5, 5, 2, 2)))
      pushViewport(dataViewport(xscale = xlim, yscale = c(-125, 125)))
      
      grid.rect()
      tmp <- pretty(xlim)
      xlabels <- tmp[tmp >= xlim[1] & tmp <= xlim[2]]
      grid.xaxis(at = xlabels)
      grid.yaxis(at = c(-100, -50, 0, 50, 100))
      
      grid.segments(top_plot$mz,
                    top_plot$intensity,
                    top_plot$mz,
                    rep(0, length(top_plot$intensity)),
                    default.units = "native",
                    gp = gpar(lwd = 0.75, col = "blue"))
      
      grid.segments(bottom_plot$mz,
                    -bottom_plot$intensity,
                    bottom_plot$mz,
                    rep(0, length(bottom_plot$intensity)),
                    default.units = "native",
                    gp = gpar(lwd = 0.75, col = "red"))
      
      grid.abline(intercept = 0, slope = 0)
      
      grid.text("intensity (%)", x = unit(-3.5, "lines"), rot = 90)
      grid.text("m/z", y = unit(-3.5, "lines"))
      
      popViewport(1)
      pushViewport(dataViewport(xscale = c(0, 20), yscale = c(-10, 10)))
      grid.text(top.label, unit(10, "native"), unit(9, "native"))
      grid.text(bottom.label, unit(10, "native"), unit(-9, "native"))
      
      popViewport(2)
      
    }
    
    p <- grid.grabExpr(headTailPlot())
    
  }
  
  if(output.list == FALSE) {
    
    return(similarity_score)
    
  }
  
  if(output.list == TRUE) {
    
    return(list(similarity.score = similarity_score, 
                alignment = alignment,
                plot = p))
    
  }
  
  # simscore <- as.vector((u %*% v)^2 / (sum(u^2) * sum(v^2)))   # cos squared
  
}
