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

**S'ha de crear un script per a generar el** ```taxon-set.conf```
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



























