#!/bin/bash

# Usage: ./count.sh $1
# $1 = alignments 

for i in $1/*.fasta; do
        file=$(basename $i)
        uce=${file%.fasta}
        sed -i "s.>.>${uce}_.g" $i      
done
