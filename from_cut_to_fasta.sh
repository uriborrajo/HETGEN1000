#!/bin/bash

for i in $1/*; do
        file=$(basename $i)
        uce_fasta=${file%.cut}
        mv $1/$file $1/$uce_fasta
done
