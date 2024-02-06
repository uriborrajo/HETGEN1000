#!/bin/bash

# Usage: ./zorro_mask.sh $1 $2
# $1 = alignments (UCEs)
# $2 = output directory

mkdir -p $2

for i in $1/*.fasta; do

        /home/intern/Desktop/apps/zorro_linux_x86_64 -sample $i > ${i}.mask 

done

mv $1/*.mask $2 
