# Upload Spatial metabo data
# many thanks to  :  Henrik Heiland and Jan Kückelhaus from theMILOlab
# https://github.com/theMILOlab/Code_Request/blob/abd61a05f634333975e8342220633faaac8929fa/MALDI.R
# 

# https://www.bioconductor.org/packages/release/bioc//vignettes/Cardinal/inst/doc/Cardinal-3-guide.html
library(Cardinal)
mywdir <- "~/spatial_thesis/data_and_analyses/spatial_metabolomics/"
file_mass <- "datametabo_katia/test3-gbm-lac-dan-1-1/test3-gbm-lac-dan-1-1.imzML"
m <- readMSIData(paste0(mywdir,file_mass))
m

### note : a list with the m/z and molecules names is missing
# peaks m/z de la base de donnees  MassBank:
# https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-Keio_Univ-KO000805
# https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-NAIST-KNA00278
# https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-Keio_Univ-KO001271
# 

# plotting spectra of one pixel in coordinates y=700 z=700 : 
plot(m, pixel=c(700,700))  


# plotting metabolites intensities by their m/z peak: 
image(m, mz=73,
      colorkey = FALSE) # glucose ? sans colorkey

image(m, mz=73) # glucose ?

image(m, mz=73) # glucose ?

image(m, mz=59.22) # glucose ?

image(m, mz=59.22,
      colorkey = FALSE) # glucose ?

image(m, mz=129) # glutamate ?

image(m, mz=146.045)  # glutamate ?

image(m, mz=89.0) # lactate ?

# pdf(file=paste0(mywdir, "results/xx.pdf"), width=90, height=40)
# image(m, mz=45)
# dev.off()    # -> pdf is low quality


plot(mz(m),featureData(m),type='pch')
#Error in xy.coords(x, y, xlabel, ylabel, log) : 
#  'x' and 'y' lengths differ