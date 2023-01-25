

#######################################################
################ visualize spatial metabolomics
# use cardinal to visualize spatial metabolomics : 
# and use also https://paoloinglese.github.io/tut_maldiquant_msi/
library(Cardinal)
library(MALDIquant)
library(MALDIquantForeign)
library(irlba)
library(viridis)


f = "visium_coauthorsBMC/doi_10.5061_dryad.h70rxwdmj__v11/MALDI_1/MALDI_1/raw/20201029scilslab_ncfr_glia_combine_root_mean_square"
print(file.exists(paste0(f,".imzML")))
smetabo = Cardinal::readImzML(file=paste0(f, ".imzML"), name=f) # demands ibd file! 



# note that authors merged the peaks and metadata in the same imzML

folde="visium_coauthorsBMC/doi_10.5061_dryad.h70rxwdmj__v11/MALDI_1/MALDI_1/raw/"
datamz = readImzML(folder=folde, name="20201029scilslab_ncfr_glia_combine_root_mean_square")


pek = MALDIquantForeign::importImzMl(path=paste0(rawmef, ".imzML"),
                                     centroided=TRUE, verbose=FALSE)
