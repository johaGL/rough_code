# SPATA: install
# https://themilolab.github.io/SPATA/articles/spata-installation.html
# $ conda activate scst
# $ rstudio  # we are here

# devtools::install_github(repo = "kueckelj/confuns") # ok
# devtools::install_github(repo = "theMILOlab/SPATA") # ERROR
#ERROR: dependencies ‘concaveman’, ‘ggalt’, ‘magick’ are not available for package ‘SPATA’
#* removing ‘/home/johanna/miniconda3/envs/scst/lib/R/library/SPATA’

# solution: ?
# conda install -c conda-forge r-concaveman
# conda install -c conda-forge r-ggalt
# conda install -c conda-forge r-magick
# conda install -c bioconda bioconductor-singlecellexperiment
devtools::install_github(repo = "theMILOlab/SPATA") # this time ok :)







