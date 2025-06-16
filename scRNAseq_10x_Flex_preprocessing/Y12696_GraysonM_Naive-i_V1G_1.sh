#!/usr/bin/bash
#SBATCH --account PAS2505
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64GB


date

CellRanger=/fs/ess/PCON0022/tools/cellranger-7.1.0/cellranger
cd /fs/ess/PAS2505/230801_Grayson_GSL-RH-3496/alignment/alignment_epi_data

${CellRanger} multi --id=Y12696_GraysonM_Naive-i_V1G_1 --csv=Y12696_GraysonM_Naive-i_V1G_1.csv --localcores=8 --localmem=64

date
 
