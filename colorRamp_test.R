library(ggplot2)
library(RColorBrewer)

# how to map colors yourself
# input values
lfcs = c(1.5, 2.6, -2, -1.6, -0.8, 2.4, 1.5, 0, 0.3, -1.1,-1.6, -0.2)
df <- data.frame(x=seq(length(lfcs)), y = lfcs)

# function mymapcol
mymappedcolors <- function(inputvals, collow, colmiddle, colhigh){
  elems = length(inputvals)
  if (elems < 100){
    mycolorstomap = colorRampPalette(c(collow, colmiddle, colhigh))(100)
  }else{
    mycolorstomap = colorRampPalette(c(collow, colmiddle, colhigh))(elems+1)
  }
  myabsv = max(abs(inputvals))
  
  minmax = c(-myabsv, myabsv)
  # distribute colors intervals:
  n = length(mycolorstomap)
  colvavec = seq(from=minmax[1],to=minmax[2],length.out=n)
  colvavec
  names(colvavec) = mycolorstomap
  
  storage = c()
  for (i in 1:length(inputvals)){
    if (inputvals[i] < 0){
      co = names(which.min(abs(abs(colvavec[colvavec<0])+inputvals[i])))
    }else if (inputvals[i] > 0){
      co = names(which.min(abs(colvavec[colvavec>0]-inputvals[i])))
    }else if (inputvals[i] == 0){
      co = names(colvavec[n/2]) # the center
    }
    storage = c(storage, co)
  }
  return(storage)
}

ay = rep(c("blue","yellow"),5)
ay
storage = mymappedcolors(df$y,"darkgreen","white","orange")
ggplot(df, aes(x, y)) + geom_point(size=5, color=storage)


