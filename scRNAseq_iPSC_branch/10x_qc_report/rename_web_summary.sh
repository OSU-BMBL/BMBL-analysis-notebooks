#/usr/bin/bash

wd=/mnt/c/Users/flyku/Documents/GitHub/NOTCH1-scRNAseq/10x_qc_report/
cd $wd

files=$(find . -type f -name "web_summary.html" -printf "/%P\n") 

for FILE in $files; do
  DIR=$(dirname "$FILE")
  BASE=$(basename -- "$FILE")
  EXTENSION="${BASE##*.}"	
  echo $DIR
  cp ."$FILE" $wd"$DIR"."$EXTENSION"
done 

