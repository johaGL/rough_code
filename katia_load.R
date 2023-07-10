# Upload Spatial metabo data
# many thanks to  :  Henrik Heiland and Jan Kückelhaus from theMILOlab
# https://github.com/theMILOlab/Code_Request/blob/abd61a05f634333975e8342220633faaac8929fa/MALDI.R
# 

# https://www.bioconductor.org/packages/release/bioc//vignettes/Cardinal/inst/doc/Cardinal-3-guide.html
library(Cardinal)

file_mass <- "~/datametabo_katia/test3-gbm-lac-dan-1-1/test3-gbm-lac-dan-1-1.imzML"
m <- readMSIData(file_mass)
m
plot(m, pixel=c(700,700))

image(m, mz=70) # glucose ?
image(m, mz=215) # glucose ?
image(m, mz=45) # lactate ?
image(m, mz=146) # glutamate ?

plot(mz(m),featureData(m),type='pch')
#Error in xy.coords(x, y, xlabel, ylabel, log) : 
#  'x' and 'y' lengths differ