library(Seurat)
library(devtools)
# notes on SPATA2 installation:
# install SPATA2: https://github.com/theMILOlab/SPATA2/issues/16
# still failing: devtools::install_github(repo = "theMILOlab/SPATA2")
# When using M1 systems (my case) clone it and do as said in : 
# https://github.com/theMILOlab/SPATA2/issues/41
devtools::install("~/SPATA2")
# selected patients : 248 (patient 1) , 259 (patient 2) 275 (patient3),
# and the patients 243, 334, 313, 255

p488_spato = SPATA2::transformSeuratToSpata()
