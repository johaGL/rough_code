
env_lib <- "miniconda3/envs/r_env_tflex/lib/R/library/"
.libPaths(c(env_lib, .libPaths()))
library(ggplot2)

set.seed(6)
dat <- data.frame(x = rnorm(100), y = rnorm(100), z=rnorm(100))
ggplot(dat, aes(x,y,color=z)) + geom_point() +
  scale_color_gradientn(colours=c("whitesmoke", "khaki3", "firebrick"), limits=c(-3,3),
                        breaks=c(-3/9*8,-3/9*4,0,3/9*4,3/9*8), labels=c(-2.4,-1.2,0,1.2,2.4), na.value = "black",
                        guide=guide_colorbar(nbin=10, raster=F, barwidth=20, frame.colour=c("black"),
                                             frame.linewidth=1, ticks.colour="black",  direction="horizontal")) +
  theme(legend.position = "bottom") 
