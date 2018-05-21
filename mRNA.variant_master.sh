#!/bin/bash
# @author Slim Fourati (sxf279@case.edu)
# @version 0.1

# read input arguments
email="sxf279@case.edu"
genome=Mmul_8
acceptedGenome=("GRCh38" "Mmul_8")

while getopts :d:e:g: option
do
    case "${option}" in
	d) fastqDir=$OPTARG;;
	e) email=$OPTARG;;
	g) genome=$OPTARG
	    if [[ ! "${acceptedGenome[@]}" =~ "$genome" ]]
	    then
		echo "Invalid -g argument: choose between ${acceptedGenome[@]}"
		exit 1
	    fi;;
	\?) echo "Invalid option: -$OPTARG"
	    exit 1;;
	:)
	    echo "Option -$OPTARG requires an argument."
	    exit 1;;
    esac
done

# test that directory contains seq files
nfiles=$(find $fastqDir -name "*_1.fq.gz" | wc -l)
if [ $nfiles -lt 1 ]
then
    echo "error...empty input directory"
    exit 1
fi

# initialize directories
# remove trailing back slash 
dataDir=$(echo $fastqDir | sed -r 's|/$||g')
dataDir=$(echo $dataDir | sed -r 's|/[^/]+$||g')

# make directories for every files and move file in directory
for i in `seq 1 $nfiles`
do
    mkdir -p $dataDir/raw$i
    find $fastqDir -name "*_1.fq.gz" | \
	head -n 1 | \
	sed -r "s/_1.fq.gz/_2.fq.gz/g" | \
	xargs -i mv "{}" "$dataDir/raw$i"
    find $fastqDir -name "*_1.fq.gz" | \
        head -n 1 | \
	xargs -i mv "{}" "$dataDir/raw$i"
done

# launch genome indexing
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    genomeGenerate_slurm.sh
cmd="sbatch genomeGenerate_slurm.sh -g $genome"
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')

# launch preprocessing slurm script
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    mRNA.variant_slurm.sh
sed -ri "s|^#SBATCH --array=1-.+$|#SBATCH --array=1-${nfiles}|g" \
    mRNA.variant_slurm.sh
sed -ri "s|^#SBATCH --depend=afterok:.+$|#SBATCH --depend=afterok:${slurmid}|g" \
    mRNA.variant_slurm.sh
cmd="sbatch mRNA.variant_slurm.sh -d $dataDir -g $genome"
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')

# launch vcf combine script
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    combineVCF_slurm.sh
sed -ri "s|^#SBATCH --depend=afterok:.+$|#SBATCH --depend=afterok:${slurmid}|g" \
    combineVCF_slurm.sh
sbatch combineVCF_slurm.sh -d $dataDir -g $genome
