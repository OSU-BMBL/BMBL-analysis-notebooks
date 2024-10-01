#!/usr/bin/bash
#SBATCH --job-name fetch_sra_with_sratoolkit
#SBATCH --account PCON0100
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=128GB

date
module load python/3.9-2022.05
module load sratoolkit/2.11.2

cd $HOME
python ./fetch-with-sratoolkit.py -o ./ -t ./ ./sra_ids.json
date
