#!/bin/bash


# Script created by Oriol Borrajo on 20 November 2023
# https://github.com/uriborrajo/HETGEN1000/


# ./fastp.sh {input PATH} {output PATH}
# e.g. ./fastp.sh ~/Desktop/Oriol/fastq ~/Desktop/Oriol/clean-fastq
conda activate phyluce-1.7.2
mkdir -p $2

for i in $1/*; do
    if [ -d $i ]; then
       ssp=$(basename $i)
       echo "PROCESSING: $ssp"
		
       r1=$(find $i -type f -name *_R1.fastq.gz)
       r2=$(find $i -type f -name *_R2.fastq.gz)
       
       if [ -n r1 ] && [ -n r2 ]; then
          o_spp=$2/$spp
	  mkdir -p $o_spp
          fastp -i $r1 -I $r2 -o $o_spp/$spp-READ1.fastq -O $o_spp/$spp-READ2.fastq --trim_poly_g --correction --overrepresentation_analysis --html --thread 10 \
       fi
    fi
done
