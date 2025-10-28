#!/bin/bash

## Set up input file
variants="input/input.tsv"
window=10000

while read line; do
  KF_id=$(echo "$line" | cut -f 1)
  chr=$(echo "$line" | cut -f 3)
  coord_pos=$(echo "$line" | cut -f 4)
  gene=$(echo "$line" | cut -f 11)
  
  prefix="$KF_id-$gene-$chr-$coord_pos"
  echo "Processing $prefix"
  
  ## TODO: Get window from splice event
  coordinates=$chr":"$(($coord_pos - $window))"-"$(($coord_pos + $window))
    
  ## get file id
  crams=$(grep "$KF_id" input/manifest.tsv | grep "Aligned.out.sorted.cram" | grep -v "crai" | cut -f2)
    
  ## loop through each CRAM per patient
  ## TODO: Make select from BS_ID an option?
  for cram in $crams; do
    cram_path="../data/cavatica/projects/sicklera/pbta-and-normal-crams/$cram"
    
    # prefix: derive from file name
    prefix=$(basename "$cram" .Aligned.out.sorted.cram)
    
    echo "Converting $cram_path"
    bam_path="results/bams/${prefix}-${KF_id}-${gene}-${coordinates}.bam"
    # input_path="variants/${prefix}-${KF_id}-${gene}-${coordinates}.tsv"
      
    samtools view \
      -T ../data/hg38.fa \
      -b \
      "$cram_path" \
      "$coordinates" \
      -o "$bam_path"
    
    samtools index "$bam_path"
    
  done
done < <(tail -n +2 $variants)