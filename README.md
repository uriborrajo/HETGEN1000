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

input="/home/intern/Desktop/Oriol/UCE_clean_reads_mollusca"
output="/home/intern/Desktop/Oriol/cdhitdup"

mkdir -p "/home/intern/Desktop/Oriol/cdhitdup"

for species_directory in "$input"/*; do
  if [ -d "$species_directory" ]; then
    species=$(basename "$species_directory")
    echo "===== Processing $species ====="

    READ1=$(find "$species_directory" -type f -name "*-READ1.fastq")
    READ2=$(find "$species_directory" -type f -name "*-READ2.fastq")

    if [ -n "$READ1" ] && [ -n "$READ2" ]; then
      output_species="$output/$species_directory"
      mkdir -p "$output_species"  
      gunzip -c -k "$READ1" "$READ2" | cd-hit-dup -u 30 -m false -i "$READ1" -i2 "$READ2" \
        -o "$output_directory/${species}-READ1.fastq" -o2 "$output_dir/${species}-READ2.fastq"

      gzip -k "$output_species/${species}-READ1.fastq"
      gzip -k "$output_species/${species}-READ2.fastq"
      
      echo "===== Duplicates removed of $species ====="
    else
      echo "No "
    fi
  fi
done

for species_dir in "$secondary_directory"/*; do
  if [ -d "$species_dir" ]; then
    find "$species_dir" -type f ! -name "*.gz" -exec rm -f {} \;
  fi
done
```
```
#!/bin/bash

mkdir -p "/home/intern/Desktop/Oriol/cdhitdup"

directorio_principal="/home/intern/Desktop/Oriol/UCE_clean_reads_mollusca"
directorio_secundario="/home/intern/Desktop/Oriol/cdhitdup"

for especie_dir in "$directorio_principal"/split-adapter-quality-trimmed/*; do
  if [ -d "$especie_dir" ]; then
    especie=$(basename "$especie_dir")
    echo "Procesando especie: $especie"

    read1=$(find "$especie_dir" -type f -name "*-READ1.fastq.gz")
    read2=$(find "$especie_dir" -type f -name "*-READ2.fastq.gz")

    if [ -n "$read1" ] && [ -n "$read2" ]; then
      gunzip -c -k "$read1" "$read2"

      read1_1=$(find "$especie_dir" -type f -name "*-READ1.fastq")
      read2_1=$(find "$especie_dir" -type f -name "*-READ2.fastq")

      if [ -n "$read1_1" ] && [ -n "$read2_1" ]; then
        output_dir="$directorio_secundario/$especie"
        mkdir -p "$output_dir"
        cd-hit-dup -u 30 -m false -i "$read1_1" -i2 "$read2_1" -o "$output_dir/${especie}-READ1.fastq" -o2 "$output_dir/${especie}-READ2.fastq"

        gzip -k "$output_dir/${especie}-READ1.fastq"
        gzip -k "$output_dir/${especie}-READ2.fastq"

        echo "Duplicados eliminados para $especie"
      else
        echo "No se encontraron archivos READ1 o READ2 para $especie"
      fi
    fi
  fi
done

for especie_dir in "$directorio_secundario"/*; do
  if [ -d "$especie_dir" ]; then
    find "$especie_dir" -type f ! -name "*.gz" -exec rm -f {} \;
  fi
done
```
```
#prueba script
#!/bin/bash

mkdir -p "/home/intern/Desktop/Oriol/cdhitdup"

directorio_principal="/home/intern/Desktop/Oriol/UCE_clean_reads_mollusca"
directorio_secundario="/home/intern/Desktop/Oriol/cdhitdup"

for especie_dir in "$directorio_principal"/*; do
  if [ -d "$especie_dir" ]; then
    especie=$(basename "$especie_dir")
    echo "Procesando especie: $especie"

    read1=$(find "$especie_dir" -type f -name "*-READ1.fastq.gz")
    read2=$(find "$especie_dir" -type f -name "*-READ2.fastq.gz")

    if [ -n "$read1" ] && [ -n "$read2" ]; then
      gunzip -c -k "$read1" "$read2"

      read1_1=$(find "$especie_dir" -type f -name "*-READ1.fastq")
      read2_1=$(find "$especie_dir" -type f -name "*-READ2.fastq")

      if [ -n "$read1_1" ] && [ -n "$read2_1" ]; then
        output_dir="$directorio_secundario/$especie"
        mkdir -p "$output_dir"
        cd-hit-dup -u 30 -m false -i "$read1_1" -i2 "$read2_1" -o "$output_dir/${especie}-READ1.fastq" -o2 "$output_dir/${especie}-READ2.fastq"

        gzip -k "$output_dir/${especie}-READ1.fastq"
        gzip -k "$output_dir/${especie}-READ2.fastq"

        echo "Duplicados eliminados para $especie"
      else
        echo "No se encontraron archivos READ1 o READ2 para $especie"
      fi
    fi
  fi
done

for especie_dir in "$directorio_secundario"/*; do
  if [ -d "$especie_dir" ]; then
    find "$especie_dir" -type f ! -name "*.gz" -exec rm -f {} \;
  fi
done
```


### 6. SPADES
```
# script assembly.conf.pl para generar el documento assembly.conf
#!/usr/bin/perl

use strict;
use warnings;

my $directory = '/home/intern/Desktop/Oriol/cdhitdup';
my $conf_file = 'assembly.conf';

opendir(my $dh, $directory) or die "Cannot open directory: $!";
my @subdirectories = grep { !/^\./ && -d "$directory/$_" } readdir($dh);
closedir($dh);

open(my $fh, '>', $conf_file) or die "Cannot open file: $!";

print $fh "[samples]\n";

foreach my $subdir (@subdirectories) {
    my $subdir_path = "$directory/$subdir";
    print $fh "$subdir:$subdir_path\n";
}
close($fh);

print "Se ha generado el archivo $conf_file.\n";
```
```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --memory 20000 \
    --cores 30 \
```
```
#!/bin/bash

directorio_principal="/home/intern/Desktop/Oriol/cdhitdup/spades-assemblies"

directorio_clean="${directorio_principal}/clean"

mkdir -p "${directorio_clean}"

for carpeta_especie in "${directorio_principal}"/*/; do
    nombre_especie=$(basename "${carpeta_especie}")
    nombre_especie_sin_spades="${nombre_especie/_spades/}"

    if [ -e "${carpeta_especie}/contigs.fasta" ]; then
        mv "${carpeta_especie}/contigs.fasta" "${directorio_clean}/${nombre_especie_sin_spades}.contigs.fasta"
        echo "Renombrado y movido ${nombre_especie_sin_spades}.contigs.fasta"
    else
        echo "No se encontró contigs.fasta en ${nombre_especie}"
    fi
done

echo "Proceso completado"
```

-p, --parents     no error if existing, make parent directories as needed


### 7. CD-HIT
```
#!/bin/bash

directorio_principal="/home/intern/Desktop/Oriol/spades-assemblies"
directorio_secundario="/home/intern/Desktop/Oriol/log-cd-hit"
directorio_cd_hit="/home/intern/Desktop/Oriol/cd-hit"

mkdir -p "$directorio_secundario"
mkdir -p "$directorio_cd_hit"

for especie_contig in "$directorio_principal"/*; do
  if [ -f "$especie_contig" ]; then
    especie=$(basename "$especie_contig")
    echo "Procesando especie: $especie"
    
    input_file="$directorio_principal/$especie"
    output_dir="$directorio_secundario"

    cd-hit -i "$input_file" -o "$output_dir/$especie" -c 0.9 -aS 0.8 -M 50000 -T 30

    mv "$output_dir/$especie" "$directorio_cd_hit/"
   
  fi
done
```
```
http://www.bioinformatics.org/cd-hit/cd-hit-user-guide

-i  input input filename in fasta format, required
-o  output filename, required
-c  sequence identity threshold, default 0.9
this is the default cd-hit's "global sequence identity"
calculated as:
    number of identical amino acids in alignment
    divided by the full length of the shorter sequence
-G  use global sequence identity, default 1
    if set to 0, then use local sequence identity, calculated as :
    number of identical amino acids in alignment
    divided by the length of the alignment
    NOTE!!! don't use -G 0 unless you use alignment coverage controls
    see options -aL, -AL, -aS, -AS
-b  band_width of alignment, default 20
-M  max available memory (Mbyte), default 400
-n  word_length, default 5, see user's guide for choosing it
-l  length of throw_away_sequences, default 10
-t  tolerance for redundance, default 2
-d  length of description in .clstr file, default 20
    if set to 0, it takes the fasta defline and stops at first space
-s  length difference cutoff, default 0.0
    if set to 0.9, the shorter sequences need to be
    at least 90% length of the representative of the cluster
-S  length difference cutoff in amino acid, default 999999
    if set to 60, the length difference between the shorter sequences
    and the representative of the cluster can not be bigger than 60
-aL alignment coverage for the longer sequence, default 0.0
    if set to 0.9, the alignment must covers 90% of the sequence
-AL alignment coverage control for the longer sequence, default 99999999
    if set to 60, and the length of the sequence is 400,
    then the alignment must be >= 340 (400-60) residues
-aS alignment coverage for the shorter sequence, default 0.0
    if set to 0.9, the alignment must covers 90% of the sequence
-AS alignment coverage control for the shorter sequence, default 99999999
    if set to 60, and the length of the sequence is 400,
    then the alignment must be >= 340 (400-60) residues
-B  1 or 0, default 0, by default, sequences are stored in RAM
    if set to 1, sequence are stored on hard drive
    it is recommended to use -B 1 for huge databases
-p  1 or 0, default 0
    if set to 1, print alignment overlap in .clstr file
```

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

Then, we need to decide which taxa we want in our *taxon set*. So, we create a Perl script that will generate the *taxon-set.conf* with the list of those taxa that we want for our analysis.

To create the Perl script, use the following command: ```vim taxon-set.conf.pl```

 ```
#!/usr/bin/perl

use strict;
use warnings;

my $directory = '/home/intern/Desktop/Oriol/UCE_mollusca/cd-hit';
my $conf_file = 'taxon-set.conf';

opendir(my $dh, $directory) or die "Cannot open directory: $!";
my @files = grep { !/^\./ && -f "$directory/$_" } readdir($dh);
closedir($dh);

open(my $fh, '>', $conf_file) or die "Cannot open file: $!";

print $fh "[all]\n";

foreach my $file (@files) {
    my ($name) = $file =~ /^(.*)\.contigs.fasta$/;
    if ($name) {
        print $fh "$name\n";
    }
}

close($fh);

print "Se ha generado el archivo $conf_file.\n";
```
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
    --alignments mafft-nexus-edge-trimmed-gblocks \
    --output mafft-nexus-edge-trimmed-gblocks-clean \
    --cores 35 \
    --log-path log
```
**We are using ```mafft-nexus-edge-trimmed-gblocks``` but you can also use ```mafft-nexus-internal-trimmed-gblocks```, depending on the decision you made in step 11.**

### 13. FINAL DATA MATRICES
Make sure that you are in the correct directory ```~/taxon-sets/all```
```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-trimmed-gblocks-clean \
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
    --alignments mafft-nexus-edge-trimmed-gblocks-clean-50p \
    --output mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree \
    --phylip \
    --log-path log
```
### 15. DOWNSTREAM ANALYSIS
- #### IQTree
Make sure that you are in the correct directory ```~/taxon-sets/all/mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree```
```
iqtree -st DNA -ninit 10 -bb 1500 -s mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree.phylip
-pre iqtree-GHOST-50p -m GTR+FO*H4 -rcluster 10 -mrate G,R,E
```
```
iqtree -st DNA -ninit 10 -bb 1500 -s mafft-nexus-internal-trimmed-gblocks-clean-50p-IQTree.phylip -sp mafft-nexus-internal-trimmed-gblocks-clean-50p-IQTree.charsets -pre iqtree-PART-50p -m MFP+MERGE -rcluster 10 -mrate G,R,E
```
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
while read i; do cd /home/intern/Desktop/data/UCE_clean_reads/${i}/split-adapter-quality-trimmed; python ../../../scripts/Oriol/fastqCombinePairedEnd.py ${i}-READ1.fastq ${i}-READ2.fastq separator; done < especies_sin_archivos.txt
```








































