#!/usr/bin/bash
#SBATCH --account PAS2584
#SBATCH --time=10:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --mem=64GB

date

CellRanger=/fs/ess/PCON0022/tools/cellranger-7.1.0/cellranger
cd /fs/ess/PAS2584/scRNA.analysis.Jia.Apr2024/Raw.fastq.data

${CellRanger} multi --id=Y11789_LockeL_5-Control-4-Mock_V1G_1 --csv=Y11789_LockeL_5-Control-4-Mock_V1G_1.csv --localcores=8  --localmem=64

date

# sbatch xxx.sh s