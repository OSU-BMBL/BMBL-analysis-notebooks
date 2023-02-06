#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=00:30:00
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --mem=32GB

#### Set working directory and reference below: (Remember to remove the ending foreslash '/')
wd=/fs/scratch/PCON0022/workingdirectory
###############################################
### Choose reference genome
tools=/fs/ess/PCON0022/tools
gtf=$tools/genome/Homo_sapiens.GRCh38.99.gtf
#gtf=$tools/genome/Mus_musculus.GRCm38.99.gtf

###############################################

bam_dir=$wd/alignment_out

module load samtools
module load hisat2

cd $bamdir
ls
bamfiles="$(find $bam_dir -maxdepth 2 -name "*.sorted.bam" -print)"

$tools/subread/bin/featureCounts -T 8 -g gene_name --primary -a $gtf -o $wd/result/out.txt $bamfiles

sed '1d' $wd/result/out.txt | cut -f2,3,4,5,6 --complement > $wd/result/out2.txt
# remove all large alignment files when finished
#rm $bam_dir/*