#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=combineVCF

# standard output file
#SBATCH --output=combineVCF.log

# number of nodes and processors, memory required
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# time requirements
#SBATCH --time=96:00:00

# dependencies
#SBATCH --depend=afterok:$SLURM_JOB_ID

# launch executable script
while getopts d:g: option
do
    case "$option" in
        d) dirData=$OPTARG;;
        g) genome=$OPTARG;;
    esac
done

bash combineVCF.sh -d $dirData -g $genome
