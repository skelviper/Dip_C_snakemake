#!/bin/bash

#usage: ./runCHARM.sh

cd ../
mkdir -p slurm_log
snakemake --cluster 'sbatch --qos=medium --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J DipC!' --jobs 188 --resources nodes=188 --rerun-incomplete -s ./Dip_C_snakemake/Dip_C.smk --keep-going

mkdir -p ./analysis
cp CHARM/stat.ipynb ./analysis/stat.ipynb
