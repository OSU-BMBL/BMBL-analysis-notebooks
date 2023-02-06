#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=02:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --mem=64GB

###############################################
#### Set working directory and reference below: (Remember to remove the ending foreslash '/')
wd=/fs/scratch/PCON0022/workingdirectory

###############################################
### Choose which reference genome to use
tools=/fs/ess/PCON0022/tools
ref_index=$tools/genome/Homo_sapiens.GRCh38.99
#ref_index=$tools/genome/Mus_musculus.GRCm38.99
###############################################

module load samtools
module load hisat2

cd $wd
echo $NAME

mkdir $wd/log
mkdir $wd/fastp_out
mkdir $wd/result
mkdir $wd/result/pre_alignment
mkdir $wd/alignment_out

# FASTQ quality control, trimming, filtering
$tools/fastp -w 16 -i $R1 -I $R2 -o $wd/fastp_out/$NAME.R1.fastq.gz -O $wd/fastp_out/$NAME.R2.fastq.gz -h $wd/result/pre_alignment/$NAME.html -j $wd/result/pre_alignment/$NAME.fastp.json

# Reads alignment to reference genome using HISAT2
hisat2 -p 16 -x $ref_index --new-summary -1 $wd/fastp_out/$NAME.R1.fastq.gz -2 $wd/fastp_out/$NAME.R2.fastq.gz -S $wd/alignment_out/$NAME.sam 

# convert SAM file to BAM file
samtools view -S -b $wd/alignment_out/$NAME.sam -@ 16 > $wd/alignment_out/$NAME.bam

# sort bam files
samtools sort -@ 16 $wd/alignment_out/$NAME.bam -o $wd/alignment_out/$NAME.sorted.bam

# when finished, submit run_quantification.sh to generate count matrix