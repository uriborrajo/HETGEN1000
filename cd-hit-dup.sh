#!/bin/bash

# Usage: ./cd-hit-dup.sh $1 $2
# $1 = clean-fastq
# $2 = cdhitdup

#### a. Specify input and output path
input="$1"
output="$2"
mkdir -p "$output"

#### b. Loop
for species_directory in "$input"/*; do
  if [ -d "$species_directory" ]; then 
    species=$(basename "$species_directory")
    echo "PROCESSING: $species"

    read1=$(find "$species_directory" -type f -name "*-READ1.fastq") 
    read2=$(find "$species_directory" -type f -name "*-READ2.fastq") 

    if [ -n "$read1" ] && [ -n "$read2" ]; then
      output_species="$output/$species_directory" 
      mkdir -p "$output_species" #make dir if parents don't exist
      cd-hit-dup -u 30 -m false -i "$read1" -i2 "$read2" -o "$output_directory/${species}-READ1.fastq" -o2 "$output_dir/${species}-READ2.fastq" \
#### c. Compressing output
      gzip -k "$output_species/${species}-READ1.fastq"
      gzip -k "$output_species/${species}-READ2.fastq"

      echo "Duplicates removed from $species"
    else
      echo "[ERROR] no READ1 & READ2: $species"
    fi
  fi
done
#### d. Delete unnecessary files
for species_directory in "$output"/*; do
  if [ -d "$species_directory" ]; then
    find "$species_directory" -type f ! -name "*.gz" -exec rm -fr {} \;
  fi
done
