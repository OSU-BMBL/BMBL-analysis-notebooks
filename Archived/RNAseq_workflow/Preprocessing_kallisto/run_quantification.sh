#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=00:30:00
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --mem=32GB

data_dir="data"  # change to your data path
out_dir="out" # change to your path, for saving data

ref_fa="gencode.vM25.transcripts.fa" # change to your path, transcripts files
ref_index="gencode.vM25.transcripts.index" # change to your path, best at the same directory as ref_fa

if [[ -d $out_dir ]]
then
echo "Output directory exist."
else
    mkdir -p $out_dir
fi

# build index
if [[ -e $ref_index ]]
then
    echo "Ref directory exist."
else
    kallisto index -i $ref_index $ref_fa
fi

samples=("sample1" "sample2" "sample3") # replace to your sample
for NAME in ${samples[@]}
do
  echo "Processing ${NAME}..."
  date
  if [[ -d ${out_dir}/${NAME} ]]; then
    echo "Directory exist."
  else
    echo "Not exist."
    mkdir ${out_dir}/${NAME}
  fi

  # FASTQ quality control, trimming, filtering
  $tools/fastp -w 16 -i ${data_dir}/${NAME}_R1_001.fastq.gz -I ${data_dir}/${NAME}_R2_001.fastq.gz -o ${out_dir}/${NAME}/${NAME}_R1.fastq.gz -O ${out_dir}/${NAME}/${NAME}_R2.fastq.gz -h ${out_dir}/${NAME}/$NAME.html -j ${out_dir}/${NAME}/$NAME.fastp.json

  # Reads pseudoalignment to transcripts, output abundance.tsv, abundance.h5
  kallisto quant -i $ref_index -o ${out_dir}/${NAME} --bias ${out_dir}/${NAME}/${NAME}_R1.fastq.gz ${out_dir}/${NAME}/${NAME}_R2.fastq.gz

  rm ${out_dir}/${NAME}/*.fastq.gz
done
