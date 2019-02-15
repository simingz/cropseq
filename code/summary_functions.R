source('code/qq-plot.R')
summ_pvalues <- function(pvalues){
  library(gridExtra)
  # p values distribution histogram
  plot1 <- histogram(pvalues,col='grey',type="count",xlim=c(0,1),breaks=100, main= "p value distribution")
  # pvalues qq-plot
  plot2 <- qqunif.plot(pvalues, main="p value qq-plot")
  grid.arrange(plot1,plot2, ncol=2)
}

