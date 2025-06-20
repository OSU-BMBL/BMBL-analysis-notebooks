# global variables. 
# ls > ../fastq_list.txt

mkdir log
while read FASTQ1 FASTQ2 NAME
do
   echo $NAME
   sbatch --job-name=$NAME.run --output=./log/$NAME.out --export=R1=$FASTQ1,R2=$FASTQ2,NAME=$NAME run_primary_alignment.sh
   sleep 0.1s
done < fastq_list.txt
