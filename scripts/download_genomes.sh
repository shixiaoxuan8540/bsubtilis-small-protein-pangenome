#!/usr/bin/env bash
set -euo pipefail
mkdir -p data/genomes
test -s config/accessions.txt || { echo 'Missing approved accession list'; exit 1; }
datasets download genome accession --inputfile config/accessions.txt --include genome,gff3,protein,cds,gbff --filename data/genomes/ncbi_genomes.zip
unzip -o -q data/genomes/ncbi_genomes.zip -d data/genomes
