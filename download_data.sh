#!/bin/bash

set -e
set -o pipefail

# Download gencode v39 GTF
cd refs
GENCODE39="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz"

echo "Downloading ${GENCODE39##*/}"
curl -L $GENCODE39 | gunzip -c > gencode.v39.primary_assembly.annotation.gtf

# Download GRCh38 fasta file
REFERENCE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"

echo "Downloading ${REFERENCE##*/}"
curl -L $REFERENCE | gunzip -c > GRCh38.primary_assembly.genome.fa

