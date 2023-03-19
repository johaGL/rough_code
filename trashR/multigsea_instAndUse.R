# johaGL 2022
# first instal FFTW3 (Fourier tr)
# http://micro.stanford.edu/wiki/Install_FFTW3
# $ sudo apt-get install libfftw3-dev libfftw3-doc
# to be abl to install "qqconf" (dependency inst by default)
# and then do
BiocManager::install("multiGSEA")
BiocManager::install("graphite")
BiocManager::install("org.Hs.eg.db")

#### parenthesis : 
### important !!! : https://baderlab.github.io/CBW_Pathways_2020/lectures/Pathways_2020_Module2_lab_introduction_RI.pdf  ==> Do keep all genes ... do not remove non significant ones
# do : https://enrichmentmap.readthedocs.io/en/latest/GeneSets.html ,then 
# http://download.baderlab.org/PathwayCommons/PC2/v3/, then download : 
# Pathway Commons 2 HumanCyc.GSEA.gmt.gz                             2013-10-03 17:56  6.5K  
# Pathway Commons 2 HumanCyc.BIOPAX.owl.gz                           2013-10-03 17:57  6.6M  
# NOTE: these versions are old 2013, not ok
#### end parenthesis


library("multiGSEA")
library("org.Hs.eg.db")
library("graphite")
#library(metaboliteIDmapping) # TODO: install
data( transcriptome, package = "multiGSEA")
layers <- c("transcriptome")
odata <- initOmicsDataStructure( layer = layers)
odata$transcriptome <- rankFeatures( transcriptome$logFC, transcriptome$pValue)
names (odata$transcriptome) <- transcriptome$Symbol
head(odata)
# databases <- c("hsapiens","HumanCyc") # did not work
# graphite::pathways("hsapiens", "kegg") # worked
graphite::pathwayDatabases()
# there are wikipathways, kegg, pathbank, reactome, but not humancyc :(
# example to launch : 
pathways <- getMultiOmicsFeatures(dbs = databases,
                                  layer = layers,
                                  returnTranscriptome = "SYMBOL")
                                  
# note: if ever getting pathways with UNIPROT code, use:
AnnotationDbi::select(org.Hs.eg.db, "O60669", "SYMBOL","UNIPROT") # to get the row(s) uniprot and symbol

                                  


