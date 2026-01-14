#!/usr/bin/bash
#SBATCH --job-name fetch_sra_with_sratoolkit
#SBATCH --account PCON0100
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=128GB

date
module load python/3.12
module load sratoolkit/3.0.2
# cp ./fetch-with-sratoolkit.py $HOME
# cd $HOME
python ./fetch-with-sratoolkit.py -o ./ -t ./ ./sra_ids.json
date
