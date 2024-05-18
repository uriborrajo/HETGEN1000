# HETGEN1000

## GUIDE FOR UCE ANALYSIS

[![Instagram Follow](https://img.shields.io/badge/Follow_us-blue?logo=Instagram)](https://www.instagram.com/slug_lab/) 
[![Website](https://img.shields.io/badge/Visit_our_website-brightgreen)](https://sluglab.wixsite.com/slug-team)

![Ultra Conserved Elements](https://www.ultraconserved.org/assets/img/ultraconserved-header.png)


## Requirements :

- Phyluce version 1.7.2 or greater https://phyluce.readthedocs.io/en/latest/
- Iqtree http://www.iqtree.org
- Exabayes https://cme.h-its.org/exelixis/web/software/exabayes/
- Openmpi 4.1.0-10 https://www.open-mpi.org
- Libopenmpi-dev 4.1.0-10
- Gcc version 10.5.0 https://gcc.gnu.org
- Cd-hit-dup https://sites.google.com/view/cd-hit
- Fastp https://github.com/OpenGene/fastp
  
## GETTING STARTED

### 1. PHYLUCE INSTALLATION

Phyluce is a software package initially designed to analyze data extracted from ultra-conserved elements (UCEs) in the genomes of organisms. These UCEs are highly conserved regions of the genome, and in our case, we are studying organisms within the phylum Mollusca.

In the present study, we will use this software package to analyze data from 35 mollusc species.

To use Phyluce, you need to have Miniconda installed (in our case, Miniconda 3):
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```
If you want to install another version of Miniconda, just change the name of the ```.sh``` file in the ```wget``` command (e.g., from Miniconda3 to Miniconda2).

Next, you must initialize Miniconda3:
```
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
To verify that Miniconda is installed and working properly, close and reopen your terminal and then run ```conda list``. This should produce an output.

The installation of Miniconda3 is necessary to separate the downloaded packages into different environments, avoiding a mix in the base environment. Conda allows us to create different environments, and we will install Phyluce in an environment called ```phyluce-1.7.2```, where ```1.7.2`` is the version installed.

To install Phyluce 1.7.2 and create the new environment:

```
wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.2/distrib/phyluce-1.7.2-py36-Linux-conda.yml
conda env create -n phyluce-1.7.2 --file phyluce-1.7.2-py36-Linux-conda.yml
```

Once the software package is installed, we are ready to activate and use the corresponding environment. Note that whenever we want to use the Phyluce software, we will have to activate the environment where it is installed.

To activate the environment:

```
conda activate phyluce-1.7.2
```
or
```
source activate phyluce-1.7.2
```

### 2. IQ-TREE INSTALLATION
Inference of phylogenetic trees by the maximum likelihood method is one of the tasks for which the IQ-TREE program is often used. The program offers an accurate and time-efficient stochastic algorithm, sharing the same accuracy and time effectiveness with other phylogenetic inference programs, such as RAxML and PhyML (Nguyen et al., 2015).

To install IQ-TREE, we will use the conda command, as it allows us to download this software program with a simple command:
```
conda install -c bioconda iqtree
```
### 3. EXABAYES INSTALLATION
ExaBayes is a bioinformatics program for Bayesian phylogenetic analysis. It was especially inspired by MrBayes but also applies similar approaches from BEAST. It sets up a Markov chain Monte Carlo sampling method that allows determining parameters of the evolutionary model (e.g., branch length or substitution rates) and the posterior probability of a tree or topology (Andre J. et al., 2014).

For the installation of ExaBayes, follow the installation instructions (section 3) on the official ExaBayes website by clicking on this link: [ExaBayes Manual](https://cme.h-its.org/exelixis/web/software/exabayes/manual/manual.html#sec-3-1).

It should be noted that the installation of ExaBayes did not work with versions higher than GCC version 10, and we had to downgrade it.

Here are the steps we used for the installation:

```
sudo apt install gcc-10
sudo apt install g++-10
sudo apt-get install openmpi-bin libopenmpi-dev
```
```
mkdir -p Apps/Exabayes
cd Apps/Exabayes
wget https://cme.h-its.org/exelixis/resource/download/software/exabayes-1.5.1.tar.gz
tar -xvzf exabayes-1.5.1.tar.gz
```
```
cd exabayes-1.5.1.tar.gz
./configure --enable-mpi && make
```

### 4. CD-HIT INSTALLATION
CD-HIT is a program that applies an algorithm to reduce redundancy in sequence data, thereby enhancing the performance of advanced sequence analysis methods. It is based on predicting similarity to filter out extraneous sequences, converting them to their non-redundant counterparts (Fu L. et al., 2012).

For CD-HIT installation:
```
conda install -c bioconda cd-hit
```
For CD-HIT-DUP and CD-HIT-EST installation: 
```
conda install -c bioconda cd-hit-auxtools
```
For the present study, only CD-HIT-DUP has been used. This is used before making assemblies with SPAdes.

>CD-HIT is used after making assemblies with Spades.

>CD-HIT-EST is used for protein alignments, i.e., transcriptome alignments.


### 5. ASTRAL INSTALLATION
The ASTRAL tool is a Java-based package used for estimating a phylogenetic tree from a set of unrooted gene trees. It detects the species tree that maximizes the number of induced quartet trees present in the gene trees, using a statistical approach compatible with a multispecies coalescent model (visit ASTRAL's GitHub because a new version called ASTER has been released).

To install ASTRAL, we create a folder called "Apps," which will serve as our central location for all the extensions and programs we need during the course.

Then:º

```
mkdir -p Apps/Astral
cd Apps/Astral
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip
```
It should be noted that whenever we want to use ASTRAL, we will have to specify the path where the ASTRAL directory with the executable script of the program is located. 

### 6. FASTP INSTALLATION
Fastp is a fast and efficient tool for preprocessing DNA and RNA sequences in genomics and transcriptomics studies. It was developed to perform various sequencing data filtering and cleaning tasks automatically and with high speed. It includes a series of main functionalities: low quality sequence filtering, adapter trimming, error correction and short sequence removal.

To download this tool just run the following command:
```
conda install -c bioconda fastp
```

### 3. COUNT THE READ DATA

- To count the number of reads in a sequence file for a species, Unix tools can be used.

- This command counts reads from R1, which should be equal to R2. Otherwise, we will have to equalize the total reads of the two sequences (see [UNPAIRED FASTQ FILES](https://github.com/uriborrajo/HETGEN1000/edit/main/README.md#fastq-combine-paired-end-unpaired-fastq-files)).

- The command divides the output number by 4 to get the number of sequence reads.

```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

### 4. FASTP

To use fastp, we will have to run the script fastp.sh. This script uses as input the folder with the FASTQ files. The output can be directed to a folder called clean-fastq to indicate where the cleaned sequences will be stored.

An example of running the script would be:
```
./fastp.sh ~/Desktop/Oriol/fastq ~/Desktop/Oriol/clean-fastq
```
Note that if you copy the script directly into a new shell script, you will have to give it permissions to be executed. To give it these permissions, use:
```
chmod +x fastp.sh
```

### 5. CD-HIT-DUP
Before making the Spades assemblies we have to clean possible contaminations and duplicates from the raw data.

· Run the [cd-hit-dup.sh](https://github.com/uriborrajo/HETGEN1000/blob/main/cd-hit-dup.sh) shell script as
```
./cd-hit-dup.sh raw-data cdhitdup
```

###### FASTQ COMBINE PAIRED END (UNPAIRED FASTQ FILES)

For those unpaired FASTQ files (i.e., those files where R1 and R2 do not have the same number of reads), we will have to compare and discard the read-only data. Otherwise, the cd-hit-dup program will not work and will not detect the second reading of the sequences.

**IMPORTANT**: You have to verify that the FASTQ files use exactly four lines per sequence; otherwise, this program will not recognize these sequences.

The program's output will be three files: the first two contain the reads of the matching pairs, and the third one contains the reads of the unpaired sequences.

Usage: 

```
python fastqCombinePairedEnd.py input1 input2 separator
```

- *[fastqCombinePairedEnd.py](https://github.com/uriborrajo/HETGEN1000/blob/main/fastqCombinePairedEnd.py)*

- *input1* = LEFT  fastq or fastq.gz file (R1)

- *input2* = RIGHT fastq or fastq.gz file (R2)

- *separator*

To automate the process, we need to generate a text document (.txt) specifying the species for which no cd-hit-dup output has been generated. Then, we will run the following command to enter each folder of these species and execute the fastqCombinePairedEnd.py script to match the two reads of the sequences:

```
while read i; do \
cd /home/intern/Desktop/data/UCE_clean_reads/${i}/split-adapter-quality-trimmed; \
python ../../../scripts/Oriol/fastqCombinePairedEnd.py ${i}-READ1.fastq ${i}-READ2.fastq separator; \
done < especies_sin_archivos.txt
```

### 6. SPADES
The Spades program implemented in Phyluce requires a configuration file that includes a header specifying the samples to be assembled ([spades]) and also the name of each species and its respective path to the folder where the sequences are located.

To generate this file in a simple way, first navigate to the folder where the already cleaned cd-hit-dup sequences are located and then run the following commands:

``` 
echo "[samples]" > ../assembly.conf
```

``` 
for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/"; done >> ../assembly.conf
```

If the folder where the sequences are located is in another folder inside the species name folder, for example, inside the split-adapter-quality-trimmed folder, execute this command instead of the previous one:

``` 
for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/split-adapter-quality-trimmed/"; done >> ../assembly.conf
```
>An example of the conf file:
```
[samples]
Lomanotus_barlettai:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/Lomanotus_barlettai
Tubulophilinopsis_pilsbryi_ZMBN105153:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/Tubulophilinopsis_pilsbryi_ZMBN105153
Heterodoris_sp:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/Heterodoris_sp
Madrella_ferruginosa:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/Madrella_ferruginosa
```


After creating the configuration file, initiate the assembly with SPAdes using the following command:

```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --memory 20000 \
    --cores 30 \
```
*Keep in mind that the memory and core specifications may need adjustment based on the size of the data. In this study, with 35 specimens and 2227 UCEs captured, these parameters are set accordingly.*

### 8. FINDING UCE LOCI
Once we have assembled our contigs from the raw reads, it's time to identify the contigs that correspond to UCE loci and exclude those that do not.

The complication from this point on lies in the organization of the folders. Therefore, I recommend creating a folder structure that best suits your needs to streamline the process.

To identify the contigs corresponding to UCE loci, we need to execute the following command. In our case, we've executed it within the "contigs" folder, specifying that the probe sets are located two folders back and inside the "spades-assembly" folder:

```
phyluce_assembly_match_contigs_to_probes \
    --contigs . \
    --probes ../../spades-assembly/Probeset-70nt.fasta \
    --output uce-search-results \
```
### 9. EXTRACTING UCE LOCI
Having identified the UCE loci, our next step is to select the taxa for our analysis, compile a list of these taxa, and generate a configuration file that specifies which UCE loci are present in each taxon.

To generate a list of the taxa you want to use, you can directly extract the names from an Excel file if you have one. Alternatively, if you have a folder containing the taxa you want to include, you can use the following command:

```
ls -1 .
```

We will generate a file using the ```nano``` command and name it ```taxon-set.conf```. Once the configuration file is created, we will write ```[all]``` or the name of the taxon set that we are using on the first line. Then, we will copy the previously generated list of taxa into the following rows.

>An example of the conf file:
```
[all]
Acteocina_exilis 
Amphimeniidae_sp2 
Aplysiopsis_elegans 
Arion_vulgaris
```

To maintain good organization of our folders and files, we will create a folder to store the different taxon sets that we are creating. In our case, we will create a folder for the taxon set ```all```:

```mkdir -p taxon-sets/all```

Now we can execute the following command to generate the initial list of UCE loci enriched in each taxon:

```
phyluce_assembly_get_match_counts \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxon-sets/all/all-taxa-incomplete.conf
```
This will generate an output called ```all-taxa-incomplete.conf``, but we want it in FASTA format for further analysis.

To do this, we need to run the following command, which will generate a file called ```all-taxa-incomplete.fasta```:

```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../../spades-assemblies \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```
###### EXPLODING THE MONOLITHIC FASTA FILE
If we want to analyze the statistics for each taxon, we need to separate the ```all-taxa-incomplete.fasta`` file into individual FASTA files for each taxon, each containing its respective captured UCEs.

Run the following two commands: the first one separates the file by taxa, and the second one extracts statistics.

```
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
```
```
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```
>Example of statistics output:
```
## samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
Acteocina-exilis.unaligned.fasta,1295,656310,506.8030888030888,4.657000046699546,184,1477,501.0,12
```

### 10. ALIGNING UCE LOCI

In the case of our study we decided to align the UCE loci with no trim.

To do this, we must specify in the following command the option ```--no-trim``` and in addition, we must generate a ```log``` folder with the command ```mkdir -p log``` and
change the number of ```--taxa``` to the number of taxa they have in their taxon set.

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
For masking the highly variable, ambiguous or error-prone parts of the sequences we will use two different masking programs, Gblocks and ZORRO.

- #### 10.1 GBLOCKS

Gblocks is a bioinformatics tool that enhances the quality of multiple sequence alignments in phylogenetic analysis. It is a tool which detects and removes misaligned and ambiguous portions of a DNA, RNA or protein sequence alignments. ⁤

⁤Gblocks keeps only the areas of the sequences that are well and reliably aligned, which enhances the performance of the data used. ⁤

Phyluce implements a command to directly execute Gblocks:

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
For our study we have used three types of Gblocks, which we have called Gblocks1, Gblocks2 and Gblocks3:

***GBLOCKS1***:  --b1 0.5 --b2 0.85 --b3 4 --b4 8   # very restrictive

***GBLOCKS2***: --b1 0.5 --b2 0.5 --b3 6 --b4 6   # intermediate

***GBLOCKS3***: --b1 0.5 --b2 0.5 --b3 10 --b4 4   # very conservative

#### 10.2 ZORRO

*Step 1 (Workstation):
```
./zorro_mask.sh $1 $2
```
$1 = alignments (UCEs) (e.g. mafft-nexus-internal-trimmed)
$2 = zorro_mask

*Step 2 (Harvard):
```
conda activate biopython2
cd zorro_mask
python ../zorro.py .
```
$1 = zorro_mask 

```
mkdir zorro
mv zorro_mask/*.cut zorro
```
Then, scp intern...

*Step 3 (Workstation):
```
./from_cut_to_fasta.sh zorro
```
```
phyluce_align_filter_alignments --alignments zorro --output zorro_filter --input-format fasta --min-length 1 --log-path log
```
```
phyluce_align_remove_locus_name_from_files --alignments zorro_filter --output zorro_clean --cores 12 --log-path log
```
```
phyluce_align_get_only_loci_with_min_taxa --alignments zorro_clean --taxa X --percent 0.50 --output zorro_50p --cores 35 --log-path log
```
```
phyluce_align_concatenate_alignments --alignments zorro_50p --output zorro_50p_IQTree --phylip --log-path log
```
```
iqtree2 --seqtype DNA --ninit 10 -B 1500 -s zorro_50p_IQTree/zorro_50p_IQTree.phylip --prefix zorro-iqtree-GHOST-50p -m GTR+FO*H4 -T AUTO --rcluster 10 --mrate G,R,E
```

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
### 14. COUNT UCEs PER SPECIES (of each matrix)
```
phyluce_align_convert_one_align_to_another --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p --output mafft-fastas-internal-trimmed-gblocks-clean-50p --input-format nexus --output-format fasta --cores 12 --log-path log
```
```
./add_tag.sh mafft-fastas-internal-trimmed-gblocks-clean-50p
```
```
cd mafft-fastas-internal-trimmed-gblocks-clean-50p
cat * >> all-fastas-50p
```
```
phyluce_assembly_explode_get_fastas_file --input all-fastas-50p --output exploded-fastas-50p --by-taxon
```
```
for i in exploded-fastas-50p/*.fasta; do phyluce_assembly_get_fasta_lengths --input $i --csv; done
```
### 15. PREPARING DATA FOR DOWNSTREAM ANALYSIS
Make sure that you are in the correct directory ```~/taxon-sets/all```
```
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p-IQTree \
    --phylip \
    --log-path log
```
### 16. DOWNSTREAM ANALYSIS
- #### IQTree
Make sure that you are in the correct directory ```~/taxon-sets/all/mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree```
```
#!/bin/bash

# Script created by Oriol Borrajo on 20 November 2023
# https://github.com/uriborrajo/HETGEN1000/

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
New command:
```
iqtree2 --seqtype DNA --ninit 10 -B 1500 -s mafft-nexus-internal-trimmed-gblocks1-clean-50p-IQTree.phylip --prefix iqtree-GHOST-50p -m GTR+FO*H4 -T AUTO --rcluster 10 --mrate G,R,E
```

- #### ExaBayes
```
#!/bin/bash

# Script created by Oriol Borrajo on 20 November 2023
# https://github.com/uriborrajo/HETGEN1000/

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

When run has completed, calculate these in an interactive session:
 
- Check if parameters converged (ESS should be >100, PSRF should be <1.1)
postProcParam -f ExaBayes_parameters.17taxa.0 ExaBayes_parameters.17taxa.1 -n 17params
- Calculate standard deviation of split frequencies (aka convergence, ASDSF should be <1%)
sdsf -f ExaBayes_topologies.run-0.24Ony-6 ExaBayes_topologies.run-1.24Ony-6
- Compute credible set of topologies (aka frequency of diff topologies)
credibleSet -f ExaBayes_topologies.17taxa.0 ExaBayes_topologies.17taxa.1 -n 17Cred
- Extract bipartitions (see how well supported different nodes are)
extractBips -f ExaBayes_topologies.17taxa.0 ExaBayes_topologies.17taxa.1 -n 17Bips
[Output files to care about: .bipartitions = name of node; .bipartitionStatistics = support for node (check ESS)]
- Generate consensus tree
consense –f ExaBayes_topologies.16Taxa.0 ExaBayes_topologies.16Taxa.1 –n 16ConsTree


- #### ASTRAL
```
iqtree2 -S mafft-nexus-internal-trimmed-gblocks1-clean-50p/ --prefix loci -T AUTO --seqtype DNA -m GTR+FO*H4 -B 1500
```

```
java -jar /home/intern/Desktop/apps/ASTRAL/astral.5.7.8.jar -i *.treefile -o astral_sptree.treefile
```

## REFERENCES

• Crotty, S.M., Minh, B.Q., Bean, N.G., Holland, B.R., Tuke, J., Jermiin, L.S., Haeseler, A.V., 2019. GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol. 69, 249–264. https://doi.org/10.1093/sysbio/syz051

• Wang, H.C., Minh, B.Q., Susko, S., Roger, A.J., 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. Syst. Biol., 67:216-235. https://doi.org/10.1093/sysbio/syx068

• Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., Haeseler, A.V., Jermiin, L.S., 2017. ModelFinder: Fast Model Selection for Accurate Phyloge- netic Estimates, Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285

• Faircloth, B.C., 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. https://doi.org/10.1093/bioinformatics/btv646.

• Andre J. Aberer, Kassian Kobert, Alexandros Stamatakis, ExaBayes: Massively Parallel Bayesian Tree Inference for the Whole-Genome Era, Molecular Biology and Evolution, Volume 31, Issue 10, October 2014, Pages 2553–2556, https://doi.org/10.1093/molbev/msu236

• Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, Bui Quang Minh, IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies, Molecular Biology and Evolution, Volume 32, Issue 1, January 2015, Pages 268–274, https://doi.org/10.1093/molbev/msu300

• Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics (Oxford, England), 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565










































