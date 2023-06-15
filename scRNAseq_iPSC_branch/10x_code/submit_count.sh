# global variables. 
mkdir log
while read NAME 
do
   sbatch --job-name=$NAME.run --output=./log/$NAME.out --export=NAME=$NAME cellranger_count.sh
   sleep 0.1s
done < sample.txt
