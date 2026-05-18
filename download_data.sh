#!/bin/bash

set -e
set -o pipefail

# Download gencode v39 GTF
cd refs
GENCODE39="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/gencode.v39.primary_assembly.annotation.gtf.gz"

echo "Downloading ${GENCODE39##*/}"
curl -L $GENCODE39 > gencode.v39.primary_assembly.annotation.gtf

# Download GRCh38 fasta file
REFERENCE="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/GRCh38.primary_assembly.genome.fa"

echo "Downloading ${REFERENCE##*/}"
curl -L $REFERENCE > GRCh38.primary_assembly.genome.fa

# Download fasta index
INDEX="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/GRCh38.primary_assembly.genome.fa.fai"

echo "Downloading ${INDEX##*/}"
curl -L $INDEX > GRCh38.primary_assembly.genome.fa.fai