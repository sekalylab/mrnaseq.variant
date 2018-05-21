#!/bin/bash
# @author: Slim Fourati

# load apps
module load gatk/3.8

# read arguments
genome="Mmul_8"
while getopts d:g: option
do
    case "$option" in
	d) dataDir=$OPTARG;;
	g) genome=$OPTARG;;
    esac
done

# set global variables for the script
seqDependencies="/mnt/projects/SOM_PATH_RXS745U/genome/$genome"
genomeFasta="$seqDependencies/Sequence/genome.fa"
maxProc=8

# Combine vcf filtered vcf files
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: combining vcf..."
    vcfFiles=$(find $dataDir -name "*.filtered.vcf")
    vcfFiles=$(echo $vcfFiles | sed -r 's/ / --variant /g')
    cmd="java -d64 -jar $GATK"
    cmd="$cmd -T GenotypeGVCFs"
    cmd="$cmd -R $genomeFasta"
    cmd="$cmd -nt $maxProc"
    cmd="$cmd --variant $vcfFiles"
    cmd="$cmd -o $dataDir/variant.vcf 2>$dataDir/gvcf.log"
    # echo $cmd
    eval $cmd
    echo "done"
fi
