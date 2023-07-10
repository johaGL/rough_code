# Upload Spatial metabo data
# many thanks to  :  Henrik Heiland and Jan Kückelhaus from theMILOlab
# https://github.com/theMILOlab/Code_Request/blob/abd61a05f634333975e8342220633faaac8929fa/MALDI.R
# 
library(Cardinal)

file_m <- "~/spatialmetabo_germanteam/MALDI_raw/MALDI_CC/imzML_file/20201029scilslab_ncfr_glia_combine_root_mean_square.imzML"
m_heiland <- readMSIData(file_m)
m_heiland
plot(m_heiland, pixel=c(700,700))

image(m_heiland, mz=215) # glucose 
image(m_heiland, mz=146) # glutamate
image(m_heiland, mz=110) # ? peaks inferior to 100 do not plot
# I do not find the lactate peak in the list : Metaspace_annotations....xlsx