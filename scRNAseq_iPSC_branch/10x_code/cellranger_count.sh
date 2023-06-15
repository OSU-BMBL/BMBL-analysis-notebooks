#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=16:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=28
#SBATCH --mem=256GB

wd=/fs/scratch/PAS1475/mingtao
Collaborator=Mingtao #change

CellRanger=/users/PAS1571/wangcankun100/tools/cellranger-4.0.0/cellranger
FastqFolder=/fs/project/PAS1475/cankun/mingtao/usftp21.novogene.com/raw_data
Refer=/fs/project/PAS1475/tools/refdata-gex-GRCh38-2020-A

#########################

cd $wd
FASTQ_FOLDER=`echo $NAME | sed 's/\_.*//g'`
# conduct alignment and read count:
${CellRanger} count --id=$NAME --transcriptome=${Refer} --fastqs=${FastqFolder}/$FASTQ_FOLDER --sample=$NAME --localcores=28 --localmem=256
