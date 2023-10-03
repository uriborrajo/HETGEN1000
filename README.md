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

If you want to leave Phyluce enviroment, run:
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

# Crear directorio cdhitdup para output
mkdir -p "/home/intern/Desktop/Oriol/cdhitdup"

# Directorio principal que contiene las carpetas de especies
directorio_principal="/home/intern/Desktop/Oriol/UCE_clean_reads_mollusca"
directorio_secundario="/home/intern/Desktop/Oriol/cdhitdup"

# Itera sobre las carpetas de especies
for especie_dir in "$directorio_principal"/*; do
  if [ -d "$especie_dir" ]; then
    especie=$(basename "$especie_dir")
    echo "Procesando especie: $especie"

    # Encuentra archivos READ1 y READ2 en la carpeta de especie
    read1=$(find "$especie_dir" -type f -name "*-READ1.fastq")
    read2=$(find "$especie_dir" -type f -name "*-READ2.fastq")

    # Verifica si se encontraron archivos READ1 y READ2
    if [ -n "$read1" ] && [ -n "$read2" ]; then
      # Ejecuta CD-HIT-DUP para eliminar duplicados
      output_dir="$directorio_secundario/$especie"
      mkdir -p "$output_dir"  # Crea el directorio de salida si no existe
      gunzip -c -k "$read1" "$read2" | cd-hit-dup -u 30 -m false -i "$read1" -i2 "$read2" \
        -o "$output_dir/${especie}-READ1.fastq" -o2 "$output_dir/${especie}-READ2.fastq"

      # Comprime los archivos de salida en el directorio de especie
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







