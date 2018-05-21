#!/bin/bash
# @author: Slim Fourati

# load apps
module load STAR/2.5.3a
module load picard/2.11
module load gatk/3.8
module load samtools/1.8

# read arguments
genome="Mmul_8"
pairEnd=false
while getopts d:g: option
do
    case "$option" in
	d) dirData=$OPTARG;;
	g) genome=$OPTARG;;
    esac
done

# set global variables for the script
bin="/mnt/projects/SOM_PATH_RXS745U/bin"
seqDependencies="/mnt/projects/SOM_PATH_RXS745U/genome/$genome"
genomeFasta="$seqDependencies/Sequence/genome.fa"
maxProc=8
genomeDir="$seqDependencies/ggNoOverhang"

# Path to reference genome and Index files.
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: alignment of reads..."
    sampleID=$(find $dirData -name "*_[1-2].fq.gz")
    sampleID=$(echo $sampleID | sed -r 's/_[1-2].fq.gz//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	STAR --genomeDir $genomeDir \
	    --readFilesIn ${sample}_1.fq.gz \
	    ${sample}_2.fq.gz \
	    --readFilesCommand zcat \
            --runThreadN $maxProc \
            --outSAMtype BAM Unsorted \
	    --outSAMstrandField intronMotif \
	    --twopassMode Basic \
	    --outFileNamePrefix ${sample}_star \
	    --outSAMattrRGline ID:RhCMV \
	    CN:Gale_lab \
	    LB:PairedEnd \
	    PL:Illumina \
	    PU:Unknown \
	    SM:${sample} &>/dev/null
	if [ $? != 0 ]
        then
            echo -ne "error\n  unable to aligned read in directory $sample"
            exit
        fi
	# remove file
	rm ${sample}_1.fq.gz
	rm ${sample}_2.fq.gz
	rm ${sample}_starLog.out
	rm ${sample}_starLog.progress.out
	rm ${sample}_starSJ.out.tab
	rm -r ${sample}_star_STARgenome
	rm -r ${sample}_star_STARpass1
    done
    echo "done"
fi

# Order BAM file by position
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: sorting bam file by position..."
    for sample in $(find $dirData -name "*_starAligned.out.bam")
    do
        sampleID=$(echo $sample | sed -r 's/(.?)_starAligned.out.bam$/\1/')
	$bin/sambamba-0.6.6/sambamba_v0.6.6 sort \
	    -o ${sampleID}.sorted.bam \
	    -p \
	    -t $maxProc \
	    $sample &>/dev/null
	# remove file
	rm $sample
    done
    echo -ne "done\n"
fi

# Index BAM file
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: indexing bam file..."
    for sample in $(find $dirData -name "*.sorted.bam")
    do
	$bin/sambamba-0.6.6/sambamba_v0.6.6 index \
	    -p \
	    -t $maxProc \
	    $sample &>/dev/null
    done
    echo -ne "done\n"
fi

# Mark duplicates using Picard
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: marking duplicates..."
    for sample in $(find $dirData -name "*.sorted.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).sorted.bam$/\1/')
	java -jar $PICARD MarkDuplicates \
	    I=$sample \
	    O=$sampleID.dupMarked.bam \
	    M=$sampleID.dup.metrics \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT 2>$sampleID.markDuplicates.log
	# remove files
	rm $sample
	rm $sample.bai
    done
    echo -ne "done\n"
fi

# SplitNCigarReads
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: splitting reads..."
    for sample in $(find $dirData -name "*.dupMarked.bam")
    do
        sampleID=$(echo $sample | sed -r 's/(.?).dupMarked.bam$/\1/')
	java -d64 -jar $GATK \
	    -T SplitNCigarReads \
	    -R $genomeFasta \
	    -I $sample \
	    -o $sampleID.split.bam \
	    -rf ReassignOneMappingQuality \
	    -RMQF 255 \
	    -RMQT 60 \
	    -U ALLOW_N_CIGAR_READS 2>$sampleID.splitNCigarReads.log
	# remove files
	rm $sample
	rm $sampleID.dupMarked.bai
	rm $sampleID.dup.metrics
    done
    echo -ne "done\n"
fi

# Create targets for indel realignment
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: creating targets for realignment..."
    for sample in $(find $dirData -name "*.split.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).split.bam$/\1/')
	java -d64 -jar $GATK \
            -T RealignerTargetCreator \
	    -R $genomeFasta \
	    -I $sample \
	    -o $sampleID.intervals \
	    -nt $maxProc 2>$sampleID.indel.log
	# -known $seqDependencies/Annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf
	# -known $seqDependencies/Annotation/Homo_sapiens_assembly38.known_indels.vcf
    done
    echo -ne "done\n"
fi

# Perform indel realignment
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: perform indel realignment..."
    for sample in $(find $dirData -name "*.split.bam")
    do
        sampleID=$(echo $sample | sed -r 's/(.?).split.bam$/\1/')
        java -d64 -Xmx30G -jar $GATK \
            -T IndelRealigner \
            -R $genomeFasta \
            -targetIntervals $sampleID.intervals \
	    -I $sample \
            -o $sampleID.processed.bam 2>$sampleID.indel2.log
            # -known $seqDependencies/Annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf
            # -known $seqDependencies/Annotation/Homo_sapiens_assembly38.known_indels.vcf
	# remove files
	rm $sample
	rm $sampleID.split.bai
	rm $sampleID.intervals
    done
    echo -ne "done\n"
fi

# Generate recalibration table
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: generate recalibration table..."
    for sample in $(find $dirData -name "*.processed.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).processed.bam$/\1/')
	java -d64 -jar $GATK \
	    -T BaseRecalibrator \
	    -I $sample \
	    -R $genomeFasta \
	    -knownSites $seqDependencies/Annotation/variant.vcf \
	    -o $sampleID.recal.table 2>$sampleID.genBaseRecalib.log
    done
    echo -ne "done\n"
fi

# base recalibration
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: printing recalibrated reads..."
    for sample in $(find $dirData -name "*.processed.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).processed.bam$/\1/')
	java -d64 -Xmx30G -jar $GATK \
	    -T PrintReads \
	    -R $genomeFasta \
	    -I $sample \
	    -nct $maxProc \
	    -BQSR $sampleID.recal.table \
	    -o $sampleID.recal.bam 2>$sampleID.baseRecalib.log
	# remove files
	rm $sample
	rm $sampleID.processed.bai
    done
    echo -ne "done\n"
fi

# Variant calling
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: calling variants..."
    for sample in $(find $dirData -name "*.recal.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).recal.bam$/\1/')
	java -d64 -Xmx30G \
	    -jar $GATK \
	    -T HaplotypeCaller \
	    -R $genomeFasta \
	    -I $sample \
	    -dontUseSoftClippedBases \
	    -stand_call_conf 20 \
	    -ERC GVCF \
	    -variant_index_type LINEAR \
	    -variant_index_parameter 128000 \
	    -nct $maxProc \
	    -o $sampleID.vcf 2>$sampleID.haploCaller.log
	# -stand_emit_conf 20 is obsolete
	# remove files
	# rm $sample
	# rm $sampleID.recal.bai
	# rm $sampleID.recal.table
    done
    echo "done"
fi

# filtering variant calling
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: filtering variants..."
    for sample in $(find $dirData -name "*.vcf")
    do
	sampleID=$(echo $sample | sed -r 's/(.?).vcf$/\1/')
	java -d64 -jar $GATK \
	    -T VariantFiltration \
	    -R $genomeFasta \
	    -V $sample \
	    -window 35 \
	    -cluster 3 \
	    --filterName FS \
	    --filterExpression "FS > 30.0" \
	    --filterName QD \
	    --filterExpression "QD < 2.0" \
	    -o $sampleID.filtered.vcf 2>$sampleID.variantFilter.log
	# remove files
	# rm $sample
	# rm $sampleID.vcf.idx
    done
    echo "done"
fi
