#!/usr/bin/bash
#SBATCH --account PCON0100
#SBATCH --time=06:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --mem=128GB

###############################################
#### Set working directory and reference below: (Remember to remove the foreslash '/')
wd=/fs/ess/PCON0005/cankun/mingtao/velocity
###############################################

ml python
source activate cellrank
cd $wd
module load samtools

echo $NAME
# velocyto has bug to sort BAM files, so we need to manually run samtools before velocyto
rm $NAME/outs/cellsorted_possorted_genome_bam.bam*
samtools sort -@ 16 -t CB -O BAM -o $NAME/outs/cellsorted_possorted_genome_bam.bam $NAME/outs/possorted_genome_bam.bam
#genes.gtf is copied from 10x Human reference genome folder
velocyto run10x $NAME genes.gtf
