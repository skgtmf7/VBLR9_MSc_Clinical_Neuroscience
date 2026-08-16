#!/usr/bin/bash

#$ -N mqtl_enrichment
#$ -wd /home/ /Scratch/AD/
#$ -pe smp 4
#$ -l h_rt=48:00:00
#$ -l mem_free=4G
#$ -o ADmqtl_enrichment.out
#$ -e ADmqtl_enrichment.err

# Initialize module system
source /etc/profile

# Load R
module purge
module load r/r-4.3.3_bc-3.18

echo "Job running on $(hostname)"
echo "Working directory: $(pwd)"

# Set threads to match SGE allocation
export OMP_NUM_THREADS=$NSLOTS
export OPENBLAS_NUM_THREADS=$NSLOTS
export MKL_NUM_THREADS=$NSLOTS

echo "Cores allocated: $NSLOTS"

Rscript AD_HPC.R $NSLOTS

echo "Job finished"



