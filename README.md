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

mkdir -p "/home/intern/Desktop/Oriol/cdhitdup"

directorio_principal="/home/intern/Desktop/Oriol/UCE_clean_reads_mollusca"
directorio_secundario="/home/intern/Desktop/Oriol/cdhitdup"

for especie_dir in "$directorio_principal"/*; do
  if [ -d "$especie_dir" ]; then
    especie=$(basename "$especie_dir")
    echo "Procesando especie: $especie"

    read1=$(find "$especie_dir" -type f -name "*-READ1.fastq")
    read2=$(find "$especie_dir" -type f -name "*-READ2.fastq")

    if [ -n "$read1" ] && [ -n "$read2" ]; then
      output_dir="$directorio_secundario/$especie"
      mkdir -p "$output_dir"  
      gunzip -c -k "$read1" "$read2" | cd-hit-dup -u 30 -m false -i "$read1" -i2 "$read2" \
        -o "$output_dir/${especie}-READ1.fastq" -o2 "$output_dir/${especie}-READ2.fastq"

      gzip -k "$output_dir/${especie}-READ1.fastq"
      gzip -k "$output_dir/${especie}-READ2.fastq"
      
      echo "Duplicados eliminados para $especie"
    else
      echo "No se encontraron archivos READ1 o READ2 para $especie"
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
### 8. FINDING UCE LOCI
```
phyluce_assembly_match_contigs_to_probes \
    --contigs . \
    --probes ../../spades-assembly/Probeset-70nt.fasta \
    --output uce-search-results \
```
### 9. EXTRACTING UCE LOCI
```mkdir -p taxon-sets/all```

 ```taxon-set.conf.pl``` 
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

```
phyluce_assembly_get_match_counts \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf
```

```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../../spades-assemblies \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```
### 10. EXPLODING THE MONOLITHIC FASTA FILE
```
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
```
### 11. ALIGNING UCE LOCI

When taxa are “closely” related (< 30-50 MYA, perhaps), I think that edge-trimming alignments is reasonable. When the taxa you are interested in span a wider range of divergence times (> 50 MYA), you may want to think about internal trimming.

- #### EDGE TRIMMING
Make sure that you are in the correct directory ```~/taxon-sets/all```
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





































