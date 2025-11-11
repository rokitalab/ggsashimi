#!/bin/bash

## Define default variables
kf_id_col=1   # KF patient ID column
chr_col=3     # Chromosome
pos_col=4     # Position
label_col=11  # Additional label to add to plot for identification, i.e. gene
window=10000  # Bases to plot either side of the position given

## Set up input files
while getopts i:m:k:c:p:l: opt; do
  case "${opt}" in
    i) input_file=${OPTARG};; # header is expected
    m) manifest=${OPTARG};;
    k) kf_id_col=${OPTARG};;
    c) chr_col=${OPTARG};;
    p) pos_col=${OPTARG};;
    l) label_col=${OPTARG};;
    w) window=${OPTARG};;
  esac
done

if [[ -z "$input_file" ]]; then
  echo "Error: Input file (-i) is required."
fi

if [[ -z "$manifest" ]]; then
  echo "Error: Manifest (-m) is required."
fi

## Download genome reference
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
ref_genome="../data/hg38.fa.gz"

if [ -f "$ref_genome" ]; then
  echo "Using reference genome: $ref_genome"
else
  echo "Downloading reference genome..."
  curl -L -o "$FILE" "$URL"
  
  # Verify download succeeded
  if [ $? -eq 0 ]; then
    echo "Download complete: $FILE"
  else
    echo "Download failed."
    exit 1
  fi
fi

####################################################

# Loop through variant file
while read line; do
  KF_id=$(echo "$line" | cut -f $kf_id_col)
  chr=$(echo "$line" | cut -f $chr_col)
  coord_pos=$(echo "$line" | cut -f $pos_col)
  label=$(echo "$line" | cut -f $label_col)
  
  prefix="$KF_id-$label-$chr-$coord_pos"
  echo "Processing $prefix"
  
  ## TODO: Get window from splice event
  coordinates=$chr":"$(($coord_pos - $window))"-"$(($coord_pos + $window))
    
  ## get file id
  crams=$(grep "$KF_id" $manifest | grep "Aligned.out.sorted.cram" | grep -v "crai" | cut -f2)
    
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
      -T $ref_genome \
      -b \
      "$cram_path" \
      "$coordinates" \
      -o "$bam_path"
    
    samtools index "$bam_path"
    
  done
done < <(tail -n +2 $input_file)