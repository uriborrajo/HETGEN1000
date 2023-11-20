# HETGEN1000

## GUIDE FOR UCE ANALYSIS

### 1. PHYLUCE INSTALLATION

To perform the phylogeny of the HETGEN project, we will use the Phyluce package. Therefore, the first step will be its installation via Miniconda2.
```
conda install phyluce-1.7.2
```
### 2. ACTIVATE PHYLUCE
You need to activate the environment before you can use any of the Phyluce commands: 
```
conda activate phyluce-1.7.2
```

If you want to leave Phyluce environment, run:
```
conda deactivate
```
### 3. COUNT THE READ DATA
```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

### 4. FASTP

We will create a script that we will call ```fastp.pl``` to automate the usage of **fastp** for all species.

*Script created by Zeyuan Chen*

This script is in Perl ```.pl``` format:
```
#!/usr/bin/env perl

use strict;

my @reads_1 = glob "/dss/dsshome1/lxc00/di76tuy/Oriol/cephas-project/SRA/fastq/*_R1.fastq.gz";

my $op_dir = "clean-fastq";

`mkdir $op_dir` if (! -e $op_dir);


foreach my $reads_1 (@reads_1) {

    # Akera_bullata_ZMBN83032_R1.fastq.gz

    $reads_1 =~ /(.+)\/([^\/]+)_R1.fastq.gz/;

    my $path = $1;

    my $id = $2;

    my $id_dir = $op_dir . "/" . $id;

    `mkdir $id_dir` if (! -e $id_dir);

    my $reads_2 = $path . "/" . $id . "_R2.fastq.gz";

    my $op1 = "$id_dir/$id\_clean1.fastq.gz";

    my $op2 = "$id_dir/$id\_clean2.fastq.gz";

    system("fastp -i $reads_1 -I $reads_2 -o $op1 -O $op2 --trim_poly_g --correction --overrepresentation_analysis --html --thread 10");

}
```
To run the script:
```
perl fastp.pl
```
### 5. CD-HIT-DUP
```
#!/bin/bash

#### a. Specify input and output path
input="/home/intern/Desktop/Oriol/clean-fastq/faltantes"
output="$input/cdhitdup"
mkdir -p "$output" #make dir if parents don't exist

#### b. Loop
for species_directory in "$input"/*; do
  if [ -d "$species_directory" ]; then #Check if species directory exist.
    species=$(basename "$species_directory") #Save the name of the species as $species
    echo "================================ PROCESSING: $species ================================"

    read1=$(find "$species_directory" -type f -name "*-READ1.fastq") #Find in species directory the file (-type f) named *-READ1.fastq and save it as $READ1
    read2=$(find "$species_directory" -type f -name "*-READ2.fastq") #Find in species directory the file (-type f) named *-READ2.fastq and save it as $READ2

    if [ -n "$read" ] && [ -n "$read2" ]; then
      output_species="$output/$species_directory" #Get output path of each species as output_species
      mkdir -p "$output_species" #make dir if parents don't exist
      cd-hit-dup -u 30 -m false -i "$read1" -i2 "$read2" -o "$output_directory/${species}-READ1.fastq" -o2 "$output_dir/${species}-READ2.fastq" \
#### c. Compressing output
      gzip -k "$output_species/${species}-READ1.fastq"
      gzip -k "$output_species/${species}-READ2.fastq"

      echo "======================== Duplicates removed of $species ========================"
    else
      echo "==================== [ERROR] no READ1 & READ2: $species ===================="
    fi
  fi
done
#### d. Removing unnecessary files
for species_directory in "$output"/*; do
  if [ -d "$species_directory" ]; then
    find "$species_directory" -type f ! -name "*.gz" -exec rm -fr {} \;
  fi
done
```


### 6. SPADES

Steps to crate "assembly.conf" file:
``` cd cdhitdup ```

``` echo "[samples]" > ../assembly.conf ```

``` for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/"; done >> ../assembly.conf ```

or

``` for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/split-adapter-quality-trimmed/"; done >> ../assembly.conf ```


```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --memory 20000 \
    --cores 30 \
```
```
#!/bin/bash

input="/home/intern/Desktop/Oriol/cdhitdup/spades-assemblies"
output="${input}/clean"
mkdir -p "${output}"
for species_directory in "${input}"/*/; do
    species=$(basename "${species_directory}")
    rspecies="${species/_spades/}"

    if [ -e "${species_directory}/contigs.fasta" ]; then
        mv "${species_directory}/contigs.fasta" "${output}/${rspecies}.contigs.fasta"
        echo "========= ${nombre_especie_sin_spades}.contigs.fasta renamed and moved ========="
    else
        echo "================== [ERROR] No contigs.fasta for ${species} =================="
    fi
done
```
>### WARNING!!
>
>Based on [Brant's](https://gist.github.com/brantfaircloth/e48e7e4eb9748854962863d104f94095) python script we keep those loci that match more than one contig.
>
>```
>phyluce_assembly_match_contigs_to_probes \
>    --contigs . \
>    --probes ../../spades-assembly/Probeset-70nt.fasta \
>    --output uce-search-results \
>    --keep-duplicates duplicates.txt
>```
>```
>mkdir -p duplicates
>mv duplicates.txt duplicates/
>cd duplicates/
>```
>```
>python ./phyluce_assembly_parse_duplicates_file.py \
>    --contigs ../ \
>    --duplicates-file duplicates.txt \
>    --output duplicates.fasta \
>    --exclude-cnt 2
>```
>```mv duplicates.fasta ../taxon-sets/all```
>```
>phyluce_assembly_explode_get_fastas_file \
>    --input duplicates.fasta \
>    --output exploded-fastas \
>    --by-taxon
>```
>```cd exploded-fastas```
>
>```cat *-DUPE1.unaligned.fasta >> duplicates.fasta ```
>
>```sed -i 's/_DUPE1//g' duplicates.fasta```
>
>```cat all-taxa-incomplete.fasta duplicates.fasta >> all-taxa-incomplete2.fasta```


### 8. FINDING UCE LOCI
We want to locate which CONTIGS are in a UCE loci and remove those that are not. Hence, we need the Probeset that is in the folder *../../spades-assembly/* and is called *Probeset-70nt.fasta*. 

**The probeset was provided by Dr. Juan Moles.**
```
phyluce_assembly_match_contigs_to_probes \
    --contigs . \
    --probes ../../spades-assembly/Probeset-70nt.fasta \
    --output uce-search-results \
```
### 9. EXTRACTING UCE LOCI
Now that we have located UCE loci, we need to determine which taxa we want to include in our analysis.

First, we make de directory where we will use to locate the output of the "data matrix configuration file".

```mkdir -p taxon-sets/all```

Then, we need to decide which taxa we want in our *taxon set*. Create taxon-set.conf with the taxa we want in our analysis.

Once we have created the conf file, we run the following command to generate the initial list of UCE loci we enriched in each taxon:

*In the command, we specify that the loci are two folders back because we are in the '/all' directory. However, if you are in a different directory, you'll need to change it to the correct path.*

```
phyluce_assembly_get_match_counts \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf
```
Finally, we have to use that list to extract the FASTA data for each taxon for each UCE loci:

```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../../spades-assemblies \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```

*The extracted FASTA data are in a monolithic FASTA file (all data for all organisms) named all-taxa-incomplete.fasta.*

### 10. EXPLODING THE MONOLITHIC FASTA FILE

Sometimes, we want to know the individual statistics on UCE assemblies for each taxon. We can do that, exploding the monolithic FASTA file into a file with the UCE loci we have enriched by taxon:

```
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
```
Then, run the stats on those exploded files:

```
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```
*samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb*

*Acteocina-exilis.unaligned.fasta,1295,656310,506.8030888030888,4.657000046699546,184,1477,501.0,12*


### 11. ALIGNING UCE LOCI

When taxa are "closely" related (<30-50 million years ago, perhaps), I consider edge-trimming alignments to be a reasonable approach. However, when the taxa you are interested in have a broader range of divergence times (>50 million years), you might want to consider internal-trimming.

Here, we show the commands for both methods:

- #### EDGE TRIMMING
```
phyluce_align_seqcap_align \
    --input all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 33 \
    --aligner mafft \
    --cores 30 \
    --incomplete-matrix \
    --output-format fasta \
    --log-path log
```
*Replace ```--taxa``` with the number of taxa you have in your analysis. Furthermore, we specify ```--log-path``` to indicate the path where we will save the command's log, so you need to create that directory using ```mkdir log```.*

- #### INTERNAL TRIMMING
```
phyluce_align_seqcap_align \
    --input all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 33 \
    --aligner mafft \
    --cores 35 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
```
*Replace ```--taxa``` with the number of taxa you have in your analysis. Furthermore, we specify ```--log-path``` to indicate the path where we will save the command's log, so you need to create that directory using ```mkdir log```.*

- #### GBLOCKS
We are going to trim this loci using **Gblocks**

So we can:
- Run gblocks trimming on the edge trimmed alignments
 ```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft-nexus-edge-trimmed \
        --output mafft-nexus-edge-trimmed-gblocks \
        --b1 0.5 \
        --b2 0.85 \
        --b3 4 \
        --b4 8 \
        --cores 12 \
        --log-path log
```
- Run gblocks trimming on the internal trimmed alignments
```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft-nexus-internal-trimmed \
        --output mafft-nexus-internal-trimmed-gblocks \
        --b1 0.5 \
        --b2 0.85 \
        --b3 4 \
        --b4 8 \
        --cores 12 \
        --log-path log
```
*higher level:  --b1 0.5 --b2 0.85 --b3 4 --b4 8  #Very conservative*

*mid level: --b1 0.5 --b2 0.5 --b3 6 --b4 6 #This is what I start with, and use in most pubs, higher-level*

*species level: --b1 0.5 --b2 0.5 --b3 10 --b4 4  #Use this with shallow datasets (species- and population-level)*

### 12. ALIGNMENT CLEANING
Make sure that you are in the correct directory ```~/taxon-sets/all```
```
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 35 \
    --log-path log
```
**We are using ```mafft-nexus-edge-trimmed-gblocks``` but you can also use ```mafft-nexus-internal-trimmed-gblocks```, depending on the decision you made in step 11.**

### 13. FINAL DATA MATRICES
Make sure that you are in the correct directory ```~/taxon-sets/all```
```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 33 \
    --percent 0.50 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --cores 35 \
    --log-path log
```
### 14. PREPARING DATA FOR DOWNSTREAM ANALYSIS
Make sure that you are in the correct directory ```~/taxon-sets/all```
```
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p-IQTree \
    --phylip \
    --log-path log
```
### 15. DOWNSTREAM ANALYSIS
- #### IQTree
Make sure that you are in the correct directory ```~/taxon-sets/all/mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree```
```
#!/bin/bash

## ./iqtree.sh {PATH} *.phylip iqtree-GHOST-50p
cd "$1"
echo "ENTERING: $1" 
iqtree -st DNA -ninit 10 -bb 1500 -s "$2" -pre "$3" -m GTR+FO*H4 -rcluster 10 -mrate G,R,E

## ./iqtree.sh {PATH} *.phylip *.charsets iqtree-PART-50p
# cd "$1"
# echo "ENTERING: $1" 
# iqtree -st DNA -ninit 10 -bb 1500 -s "$2" -sp "$3" -pre "$4" -m MFP+MERGE -rcluster 10 -mrate G,R,E
```
Old command:
```
iqtree -st DNA -ninit 10 -bb 1500 -s mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree.phylip -pre iqtree-GHOST-50p -m GTR+FO*H4 -rcluster 10 -mrate G,R,E
```

- #### ExaBayes
```
#!/bin/bash

## ./exabayes.sh {PATH} *.phylip config.nex
cd "$1"
echo "ENTERING: $1" 
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run1 -s 1234 -M 1 #exabayes run1
# Total walltime elapsed: 43:30:4.35 (hh:mm:ss) -M 3-
# Total CPU time elapsed: 1044:01:44.38 (hh:mm:ss) -M 3-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run2 -s 1234 -M 1 #exabayes run2
# Total walltime elapsed: 39:11:50.19 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 940:44:4.51 (hh:mm:ss) -M 1-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run3 -s 1234 -M 1 #exabayes run3
# Total walltime elapsed: 27:06:23.76 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 650:33:30.29 (hh:mm:ss) -M 1-
mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -n run4 -s 1234 -M 1 #exabayes run4
# Total walltime elapsed: 25:51:55.34 (hh:mm:ss) -M 1-
# Total CPU time elapsed: 620:46:8.13 (hh:mm:ss) -M 1-
# mpirun exabayes -np 16 -R 4 -C 4 -f "$2" -m DNA -c "$3" -n run1 -s 1234 -M 1 

## ./exabayes.sh {PATH} *.phylip config.nex aln.part
# cd "$1"
# echo "ENTERING: $1" 
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run1 -s 1234 -M 1 #exabayes run1 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run2 -s 1234 -M 1 #exabayes run2 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run3 -s 1234 -M 1 #exabayes run3 w/ partition file
# mpirun exabayes -np 4 -R 1 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run4 -s 1234 -M 1 #exabayes run4 w/ partition file
# mpirun exabayes -np 16 -R 4 -C 4 -f "$2" -m DNA -c "$3" -q "$4" -n run1 -s 1234 -M 1

```

config.nex:
```
begin run;
 numruns 1
 numgen 5e6
 diagfreq 5000
 samplingfreq 500
 printfreq 100
 burninproportion 0.25
 parsimonyStart  true
 printFreq 10
 numcoupledchains 4
end;
```
###### FASTQ COMBINE PAIRED END (UNPAIRED FASTQ FILES)


```
Unpaired fastq files? Compare and discard single read data:
https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py
Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they have been trimmed and cleaned
WARNING! This program assumes that the fastq file uses EXACTLY four lines per sequence
Three output files are generated. The first two files contain the reads of the pairs that match and the third contains the solitary reads.
Usage: python fastqCombinePairedEnd.py input1 input2 separator
input1 = LEFT  fastq or fastq.gz file (R1)
input2 = RIGHT fastq or fastq.gz file (R2)
separator = character that separates the name of the read from the part that
    describes if it goes on the left or right, usually with characters '1' or
    '2'.  The separator is often a space, but could be another character (e.g. ‘/’). A
    space is used by default. If the sequence names do not contain two parts
    and you want to use the full name info to pair your sequences, use 'None'
    (as text) for the separator. Eg: python fastqCombinePairedEnd.py input1 input2 None
```

```
while read i; do \
cd /home/intern/Desktop/data/UCE_clean_reads/${i}/split-adapter-quality-trimmed; \
python ../../../scripts/Oriol/fastqCombinePairedEnd.py ${i}-READ1.fastq ${i}-READ2.fastq separator; \
done < especies_sin_archivos.txt
```








































