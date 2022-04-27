#!/usr/bin/bash
# johaGL
# many thanks to :
# https://www.linuxquestions.org/questions/programming-9/bash-combine-arrays-and-delete-duplicates-882286/
# https://unix.stackexchange.com/questions/163970/how-can-i-grep-a-certain-text-and-display-its-line-and-the-line-after-as-well
# https://stackoverflow.com/questions/27001788/passing-bash-variable-as-a-pattern-to-awk
# see non tum sample 1 : NonTum1 : D573T10
mouse_sam="/scratch/CBIB/bdartigues/Glioblastomes_P3_new/STAR_Mapping/Curie/MOUSE_XENOME/D573T10/Aligned.out.sam"
human_sam="/scratch/CBIB/bdartigues/Glioblastomes_P3_new/STAR_Mapping/Curie/HUMAN_XENOME/D573T10/Aligned.out.sam"
fqdir="/mnt/cbib/Glioblastomes_P3_new/data/Curie/data/"
odir="TinyData/"

mousechr="chr16"
humanchr="chr21"
mousetxt=${odir}"mouse16.txt"
humantxt=${odir}"human21.txt"
cids=${odir}"combiuni16m21h.txt"

# #####################################################
# save selected chromosomes to txt files
# #####################################################

#wc -l $mouse_sam # DO NOT run, too long
#wc -l $human_sam # DO NOT run, too long
#head -n 70 ${mouse_sam}
# NOTE!!: SN:chr is in header, exclude
if [ ! -f ${mousetxt} ]; then
    echo "extracting reads that aligned to MOUSE chromosome 16..."
    #grep $mousechr ${mouse_sam} | grep -v SN:chr > $mousetxt
    echo "finished mouse chr16 extraction"
    echo "extracting reads that aligned to HUMAN chromosome 21..."
    #grep $humanchr ${human_sam} | grep -v SN:chr > $humantxt
    echo  "ended chr21 human extraction"
else
    echo "already created extracted files "$mousetxt" and "$humantxt
fi

if [ ! -f ${odir}"flags_status.txt" ];then
	echo "exploring reads flags"

	#echo "Mouse flags" >> $odir"flags_status.txt"
	#awk '{ count[$2]++ } END { for (w in count) print w, count[w] }' $mousetxt >> $odir"flags_status.txt"

	#echo "Human flags" >> $odir"flags_status.txt"
	#awk '{ count[$2]++ } END { for (w in count) print w, count[w] }' $humantxt >> $odir"flags_status.txt"
fi

# #####################################################
# save reads' ids and use them to create new fastq R1 and R2 files
# #####################################################

echo "filtering uniquely mapped with NH:i:1"
#grep NH:i:1 $mousetxt > ${odir}"uniqmouse16.txt"
#grep NH:i:1 $humantxt > ${odir}"uniqhuman21.txt"
echo "extracting list of ids"
if [ ! -f $cids ];then
	declare -a li1=($(cut -f 1 ${odir}"uniqmouse16.txt"))
	declare -a li2=($(cut -f 1 ${odir}"uniqhuman21.txt"))
	echo "create a single list with all reads ids in both files, no repeated ids"
	combined=(`for I in "${li1[@]}" "${li2[@]}"; do echo "$I" ; done | sort -du`)
	echo ${combined[5]}
	for i in ${combined[@]};do
    		echo $i >> $cids
	done
else
	echo "created already a list of selected reads' ids if file "$cids
fi

echo ""

doreadfile () {
	echo "trying to extract "$1	
	OUTFILEFIN=$2
	combined=$3
	for X in ${combined[@]}; do
		#lineyes=$(grep ${X} $1)
		grep -A 3 $X $1 >> $OUTFILEFIN
		#if [ "${lineyes}" !=  "" ];then
		#	echo ${lineyes} >> $OUTFILEFIN
			#awk -v pattern="$X" '$0 ~ ""pattern {for(i=1;i<=3;i++) {getline;print} }' $1 >> $OUTFILEFIN
		#else
		#	echo "next"
		#fi
	done
}
# note: forum said to use "^"pattern but it did not worked, 
declare -a combined=($(cat $cids)); echo "finished combined retrieval"

# preserve original files, send decompressed to my wdir:
# gunzip -c ${fqdir}"D573T10.R1.fastq.gz" > $odir"D573T10.R1.fastq"
# gunzip -c ${fqdir}"D573T10.R2.fastq.gz" > $odir"D573T10.R2.fastq"

# call function with arguments : fastqfilepath outputfile listofreadsids	
#doreadfile $odir"D573T10.R1.fastq" $odir"R1_16m21h.txt" $combined
doreadfile $odir"D573T10.R2.fastq" $odir"R2_16m21h.txt" $combined



