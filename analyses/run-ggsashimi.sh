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
  input_path="results/${KF_id}-${gene}-${coordinates}.tsv"
  
  ## loop through each CRAM per patient
  ## TODO: Make select from BS_ID an option?
  for cram in $crams; do
    cram_path="../data/cavatica/projects/sicklera/pbta-and-normal-crams/$cram"
    
    # prefix: derive from file name
    prefix=$(basename "$cram" .Aligned.out.sorted.cram)
    
    echo "Converting $cram_path"
    bam_path="results/bams/${prefix}-${KF_id}-${gene}-${coordinates}.bam"
    
    samtools view \
      -T ../data/hg38.fa \
      -b \
      "$cram_path" \
      "$coordinates" \
      -o "$bam_path"
    
    samtools index "$bam_path"
    
    # create input tsv for ggsashimi
    echo "$KF_id"$'\t'"$bam_path"$'\t'"$prefix" >> "$input_path"
  done
  
  ## Prepare controls
  ## TODO: Auto-select controls, or select based on histologies
  controls=("SRR26129063" "SRR26129064" "SRR26129066")
  control_manifest="cavatica/projects/sicklera/pbta-and-normal-crams/manifest_20250825_143150.tsv"

  control_crams=$(grep -E "$(printf '%s|' "${controls[@]}" | sed 's/|$//')" "$control_manifest" | grep -v "crai" | cut -f2)
  ## loop through each CRAM
  for cram in $control_crams; do
    cram_path="cavatica/projects/sicklera/pbta-and-normal-crams/$cram"

    prefix=$(basename "$cram" .Aligned.out.sorted.cram)
    ID=$(grep "$prefix" "$control_manifest" | cut -f5 | sort -u)
    
    echo "Processing control $cram_path"
    bam_path="results/bams/${prefix}-${KF_id}-${gene}-${coordinates}.bam"

    samtools view \
      -T ../data/hg38.fa \
      -b \
      "$cram_path" \
      "$coordinates" \
      -o "$bam_path"
    
    samtools index "variants/${bam_path}"
    
    # write control path to input tsv for ggsashimi
    # Column 3 being "control" collapses all controls to one plot
    echo "$ID"$'\t'"$bam_path"$'\t'"Control" >> "$input_path"
  done
  
  # run ggsashimi
  python3 ggsashimi.py -b "$input_path" -c "$coordinates" --shrink \
      -g ../data/gencode.v39.primary_assembly.annotation.protein_coding.gtf \
      -P input/palette.txt -C 3 -O 3 -A median_j -M 3 \
      -o "plots/${gene}-${KF_id}-${coordinates}"
      
done < <(tail -n +2 $variants)