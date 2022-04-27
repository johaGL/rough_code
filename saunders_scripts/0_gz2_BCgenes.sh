#!/bin/bash
TOPL=/scratch/CBIB/jgalvis/
dataloc=${TOPL}datapub/scRNA/saunders/
mypython=${TOPL}recast_scRNApub/saunders_scripts/0_gz2_BCgenes.py
myodir=${TOPL}recast_scRNApub/saunders/

declare -a filori=(F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz \
	F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz)
declare -a decompress=(F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt \
	F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt)
declare -a filout=(frontal posterior)


cd $dataloc
# ${#array[@]} is the number of elements in the array
for ((i = 0; i < ${#filout[@]}; ++i)); do
    # bash arrays are 0-indexed
    echo " * ${filout[$i]} * "
    gunzip -k ${filori[$i]}
    grep %%CELL ${decompress[$i]} >  ${filout[$i]}barcodes.txt
    grep %%GENE ${decompress[$i]} >  ${filout[$i]}genes.txt 
    rm ${decompress[$i]}
    
    # use python script to make a unique column txt file
    python3 $mypython $dataloc ${filout[$i]}barcodes.txt
    python3 $mypython $dataloc ${filout[$i]}genes.txt
done

mv *barcodes.txt $myodir
mv *genes.txt $myodir


#gunzip -k F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz
#"grep %%CELL F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt > "frontalbarcodes.txt
#"grep %%GENE F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt > frontalgenes.txt
#rm F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt


#declare -a anb=(1 2)
#for i in ${anb[@]}; do
#	echo ${filori[$i-1]} 
#	echo ${filout[$i-1]}barcodes.txt
#	echo ${filout[$i-1]}genes.txt
#done
