#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=4:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --mem=64GB

wd=/fs/scratch/PCON0022/yjl/liuxf/HPV_mapping_results/cellranger7.1
CellRanger=~/tools/cellranger-7.1.0-virus/cellranger
FastqFolder=/fs/scratch/PCON0022/yjl/liuxf/HPV_public_data/GSE168652_RAW/CC_fastq
Refer=/fs/scratch/PCON0022/yjl/liuxf/ref_genome/HPV_18_ref_genome

#########################
date
cd $wd

${CellRanger} count --id=CC_hpv18 --transcriptome=${Refer} --fastqs=${FastqFolder} --sample=CC --chemistry=SC3P_auto --localcores=8  --localmem=64

date 