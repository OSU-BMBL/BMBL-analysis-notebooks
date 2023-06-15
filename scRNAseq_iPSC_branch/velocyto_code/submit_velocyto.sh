#!/usr/bin/bash
while read NAME
do
   echo $NAME
   sbatch --export=NAME=$NAME run_velocyto.sh
   sleep 0.1s
done < sample_list.txt
