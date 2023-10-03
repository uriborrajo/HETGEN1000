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






