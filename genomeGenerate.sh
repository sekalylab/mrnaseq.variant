#!/bin/bash
# @author Slim Fourati (sxf279@case.edu)
# @author Aarthi Talla (axt427@case.edu)
# @version 0.3

# load apps
module load gcc/6.3.0
module load STAR/2.5.3a
module load picard/2.11
module load intel/17
module load openmpi/2.0.1
module load samtools/1.8

# read input arguments
while getopts g: option
do
    case "${option}" in
	g) genome=$OPTARG;;
    esac
done

# set global variables for the script
genomeDir="/mnt/projects/SOM_PATH_RXS745U/genome/${genome}"
genomeFasta="$genomeDir/Sequence/genome.fa"
genomeDir="$genomeDir/ggNoOverhang"
maxProc=8

# 2. Generating genomes for alignment with 'STAR'
currentData=$(date +"%Y-%m-%d %X")
echo -ne "$currentData: generating genomes for alignment..."
if [ -e "$genomeDir/SA" ]
then
    rm -r $genomeDir
fi
mkdir -p $genomeDir
STAR \
    --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles $genomeFasta \
    --genomeSAsparseD 2 \
    --runThreadN $maxProc &>/dev/null
if [ $? != 0 ]
then
    echo -ne "error\n  unable to index genome"
    exit
fi
currentDir=$(pwd)
if [ -e "$currentDir/Log.out" ]
then
    rm $currentDir/Log.out
fi
# change permission of index directory
chmod g+w $genomeDir
echo "done"

# create reference genome dictionary
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: creating reference genome dictionary..."
    genomeDict=$(echo $genomeFasta | sed -r "s/.fa$|.fasta/.dict/g")
    java -jar $PICARD CreateSequenceDictionary \
        R=$genomeFasta \
        O=$genomeDict
    echo "done"
fi

# create reference genome index
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: creating genome index..."
    samtools faidx $genomeFasta
    echo "done"
fi
