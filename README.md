# HETGEN1000

## GUIDE FOR UCE ANALYSIS

[![Instagram Follow](https://img.shields.io/badge/Follow_us-blue?logo=Instagram)](https://www.instagram.com/slug_lab/) 
[![Website](https://img.shields.io/badge/Visit_our_website-brightgreen)](https://sluglab.wixsite.com/slug-team)

![Ultra Conserved Elements](https://www.ultraconserved.org/assets/img/ultraconserved-header.png)

## INDEX
### 1. [Requirements](https://github.com/uriborrajo/HETGEN1000/tree/main#1-requirements-1)
### 2. [Getting Started](https://github.com/uriborrajo/HETGEN1000/tree/main#2-getting-started-1)
#### 2.1 [Phyluce installation](https://github.com/uriborrajo/HETGEN1000/tree/main#21-phyluce-installation-1)
#### 2.2 [IQ-TREE installation](https://github.com/uriborrajo/HETGEN1000/tree/main#22-iq-tree-installation-1)
#### 2.3 [ExaBayes installation](https://github.com/uriborrajo/HETGEN1000/tree/main#23-exabayes-installation-1)
#### 2.4 [CD-HIT installation](https://github.com/uriborrajo/HETGEN1000/tree/main#24-cd-hit-installation-1)
#### 2.5 [Astral installation](https://github.com/uriborrajo/HETGEN1000/tree/main#25-astral-installation-1)
#### 2.6 [Fastp installation](https://github.com/uriborrajo/HETGEN1000/tree/main#26-fastp-installation-1)
### 3. [Phylogenomic Analysis](https://github.com/uriborrajo/HETGEN1000/tree/main#3-phylogenomic-analysis-1)
#### 3.1 [Counting the read data](https://github.com/uriborrajo/HETGEN1000/tree/main#31-counting-the-read-data-1)
#### 3.2 [Fastp](https://github.com/uriborrajo/HETGEN1000/tree/main#32-fastp-1)
#### 3.3 [CD-HIT-DUP](https://github.com/uriborrajo/HETGEN1000/tree/main#33-cd-hit-dup-1)
#### 3.4 [SPAdes](https://github.com/uriborrajo/HETGEN1000/tree/main#34-spades-1)
#### 3.5 [Finding UCE loci](https://github.com/uriborrajo/HETGEN1000/tree/main#35-finding-uce-loci-1)
#### 3.6 [Extracting UCE loci](https://github.com/uriborrajo/HETGEN1000/tree/main#9-extracting-uce-loci)
#### 3.7 [Aligning UCE loci](https://github.com/uriborrajo/HETGEN1000/tree/main#10-aligning-uce-loci)
##### 3.7.1 [Gblocks](https://github.com/uriborrajo/HETGEN1000/tree/main#101-gblocks)
##### 3.7.2 [ZORRO](https://github.com/uriborrajo/HETGEN1000/tree/main#102-zorro)
#### 3.8 [Alignment cleaning](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#38-alignment-cleaning-1)
#### 3.9 [Final data Matrices](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#39-final-data-matrices-1)
#### 3.10 [Preparing data for Downstream Analysis](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#310-preparing-data-for-downstream-analysis-1)
### 4. [Downstram Analysis](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#4-downstream-analysis)
#### 4.1 [IQ-TREE](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#41-iqtree)
#### 4.2 [ExaBayes](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#42-exabayes-1)
#### 4.3 [Astral](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#43-astral-1)
### 5. [References](https://github.com/uriborrajo/HETGEN1000/tree/main?tab=readme-ov-file#5-references-1)

## 1. REQUIREMENTS 

- Phyluce version 1.7.2 or greater https://phyluce.readthedocs.io/en/latest/
- Iqtree http://www.iqtree.org
- Exabayes https://cme.h-its.org/exelixis/web/software/exabayes/
- Openmpi 4.1.0-10 https://www.open-mpi.org
- Libopenmpi-dev 4.1.0-10
- Gcc version 10.5.0 https://gcc.gnu.org
- Cd-hit-dup https://sites.google.com/view/cd-hit
- Fastp https://github.com/OpenGene/fastp
  
## 2. GETTING STARTED

### 2.1 PHYLUCE INSTALLATION

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

### 2.2 IQ-TREE INSTALLATION
Inference of phylogenetic trees by the maximum likelihood method is one of the tasks for which the IQ-TREE program is often used. The program offers an accurate and time-efficient stochastic algorithm, sharing the same accuracy and time effectiveness with other phylogenetic inference programs, such as RAxML and PhyML (Nguyen et al., 2015).

To install IQ-TREE, we will use the conda command, as it allows us to download this software program with a simple command:
```
conda install -c bioconda iqtree
```
### 2.3 EXABAYES INSTALLATION
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

### 2.4 CD-HIT INSTALLATION
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


### 2.5 ASTRAL INSTALLATION
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

### 2.6 FASTP INSTALLATION
Fastp is a fast and efficient tool for preprocessing DNA and RNA sequences in genomics and transcriptomics studies. It was developed to perform various sequencing data filtering and cleaning tasks automatically and with high speed. It includes a series of main functionalities: low quality sequence filtering, adapter trimming, error correction and short sequence removal.

To download this tool just run the following command:
```
conda install -c bioconda fastp
```
## 3. PHYLOGENOMIC ANALYSIS

### 3.1 COUNTING THE READ DATA

- To count the number of reads in a sequence file for a species, Unix tools can be used.

- This command counts reads from R1, which should be equal to R2. Otherwise, we will have to equalize the total reads of the two sequences (see [UNPAIRED FASTQ FILES](https://github.com/uriborrajo/HETGEN1000/edit/main/README.md#fastq-combine-paired-end-unpaired-fastq-files)).

- The command divides the output number by 4 to get the number of sequence reads.

```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

### 3.2 FASTP

To use fastp, we will have to run the script fastp.sh. This script uses as input the folder with the FASTQ files. The output can be directed to a folder called clean-fastq to indicate where the cleaned sequences will be stored.

An example of running the script would be:
```
./fastp.sh ~/Desktop/Oriol/fastq ~/Desktop/Oriol/clean-fastq
```
Note that if you copy the script directly into a new shell script, you will have to give it permissions to be executed. To give it these permissions, use:
```
chmod +x fastp.sh
```

### 3.3 CD-HIT-DUP
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

### 3.4 SPADES
The SPAdes program implemented in Phyluce requires a configuration file that includes a header specifying the samples to be assembled ([spades]) and also the name of each species and its respective path to the folder where the sequences are located.

To generate this file in a simple way, first navigate to the folder where the already cleaned cd-hit-dup sequences are located and then run the following commands:

``` 
echo "[samples]" > ../assembly.conf
```

``` 
for i in *; do echo "$i:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/$i/"; done >> ../assembly.conf
```

If the folder where the sequences are located is in another folder inside the species name folder, for example, inside the split-adapter-quality-trimmed folder, execute this command instead of the previous one:

``` 
for i in *; do echo "$i:/home/intern/Desktop/data/UCE_clean_reads_cdhitdup/$i/split-adapter-quality-trimmed/"; done >> ../assembly.conf
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

### 3.5 FINDING UCE LOCI
Once we have assembled our contigs from the raw reads, it's time to identify the contigs that correspond to UCE loci and exclude those that do not.

The complication from this point on lies in the organization of the folders. Therefore, I recommend creating a folder structure that best suits your needs to streamline the process.

To identify the contigs corresponding to UCE loci, we need to execute the following command. In our case, we've executed it within the "contigs" folder, specifying that the probe sets are located two folders back and inside the "spades-assembly" folder:

```
phyluce_assembly_match_contigs_to_probes \
    --contigs . \
    --probes ../../spades-assembly/Probeset-70nt.fasta \
    --output uce-search-results \
```
### 3.6 EXTRACTING UCE LOCI
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

### 3.7 ALIGNING UCE LOCI

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

- #### 3.7.1 GBLOCKS

Gblocks is a bioinformatics tool that enhances the quality of multiple sequence alignments in phylogenetic analysis. It detects and removes misaligned and ambiguous portions of DNA, RNA, or protein sequence alignments.

Gblocks retains only the regions of the sequences that are well and reliably aligned, improving the performance and accuracy of the data used in phylogenetic studies.

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
For our study, we have used three configurations of Gblocks, which we have called Gblocks1, Gblocks2, and Gblocks3:

***GBLOCKS 1***:  --b1 0.5 --b2 0.85 --b3 4 --b4 8   # very restrictive

***GBLOCKS 2***: --b1 0.5 --b2 0.5 --b3 6 --b4 6   # intermediate

***GBLOCKS 3***: --b1 0.5 --b2 0.5 --b3 10 --b4 4   # very conservative

- #### 3.7.2 ZORRO
Zorro is a bioinformatics tool used in phylogenomics to evaluate and improve the quality of multiple sequence alignments. Unlike Gblocks, which removes misaligned regions based on predefined criteria, Zorro provides a confidence score for each position in the alignments. This score indicates the reliability of that particular position. Positions with confidence values that are too low are excluded from the alignments by Zorro, resulting in more refined alignments and, consequently, more accurate phylogenetic analyses.

To perform masking with Zorro, we will use the shell script [zorro_mask.sh](https://github.com/uriborrajo/HETGEN1000/blob/main/zorro_mask.sh), in which we need to specify the folder containing our alignments and the name of the output folder. In our case, the output folder is named ```zorro_mask```, as indicated in the following command:
```
./zorro_mask.sh $1 $2
```
***$1*** = alignments (UCEs) (e.g. mafft-nexus-internal-trimmed)

***$2*** = zorro_mask

Next we will use the Python script [zorro.py](https://github.com/uriborrajo/HETGEN1000/blob/main/zorro.py). For the process to work, we need to place this script in the general directory and run it inside the ```zorro_mask``` folder. Below you can see the commands:

```
cd zorro_mask
python ../zorro.py .
```
***$1*** = zorro_mask 

*The ```zorro_mask``` folder must contain both the ```.fasta``` and ```.fasta.mask``` files; otherwise the Python script will not work.*

*Here is an example of the files that should be in the folder:*

```
-rw-rw-r--  1 intern intern  51751 abr 17 18:42 uce-998.fasta
-rw-rw-r--  1 intern intern  25137 abr  5 00:32 uce-998.fasta.mask
```
Once the process is finished, we will create a folder named ```zorro``` where we will place the ```.cut``` files.
```
mkdir zorro
mv zorro_mask/*.cut zorro
```
Now we have a folder with only the ```.cut``` files, but for the following phyluce commands to work, we must transform these files to the *FASTA* format. To do this, we run the script [from_cut_to_fasta.sh](https://github.com/uriborrajo/HETGEN1000/blob/main/from_cut_to_fasta.sh):

```
./from_cut_to_fasta.sh zorro
```

Finally, we need to run the following commands to prepare the data for downstream analysis:

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

### 3.8 ALIGNMENT CLEANING
For downstream analysis, we need to remove the locus name from the files. Currently, each file contains both the taxon name and the locus name. By performing this step, we will simplify the file names to include only the taxon name.

To remove the locus name, use the following command:
```
phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 35 \
    --log-path log
```
### 3.9 FINAL DATA MATRICES
At this point, we are interested in minimizing noise in the data and increasing confidence in phylogenetic inferences by eliminating loci or genes that may be less informative or more susceptible to error, as well as eliminating missing data.

To achieve this, we will generate an occupancy matrix with a 50% threshold. A 50% occupancy matrix means that each column of the matrix (i.e., each locus) must be present in at least 50% of the taxa to be included. This helps in reducing the noise and improving the quality of the data for phylogenetic analysis.

Use the following command to perform this process:

```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 33 \
    --percent 0.50 \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --cores 35 \
    --log-path log
```
###### COUNT UCEs FOR EACH MATRIX
To count the UCEs retained in each species in the array, we will first need to convert the alignments to FASTA format.
```
phyluce_align_convert_one_align_to_another --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p --output mafft-fastas-internal-trimmed-gblocks-clean-50p --input-format nexus --output-format fasta --cores 12 --log-path log
```
Next, we will tag each file with the name of the corresponding UCE using the [add_tag.sh](https://github.com/uriborrajo/HETGEN1000/blob/main/add_tag.sh) script.
```
./add_tag.sh mafft-fastas-internal-trimmed-gblocks-clean-50p
```
Once all the UCEs are labeled, we need to concatenate all the files into a single monolithic file, similar to the process done previously.

Use the following commands:
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
### 3.10 PREPARING DATA FOR DOWNSTREAM ANALYSIS
To ensure that IQ-TREE, ExaBayes, and ASTRAL can recognize our files, we need to concatenate all our UCEs from the matrix into a single file in PHYLIP format.

To do this, run:
```
phyluce_align_concatenate_alignments \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-50p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-50p-IQTree \
    --phylip \
    --log-path log
```
## 4. DOWNSTREAM ANALYSIS
### 4.1 IQ-TREE
To infer the Maximum Likelihood tree, we will use IQ-TREE. For our data, we will use the GHOST model (GTR+FO*H4), but you can change it to the model that best fits your data by changing the ```-m``` (model) parameter.

Here is the command used with the new IQ-TREE (v. 2):

```
iqtree2 --seqtype DNA --ninit 10 -B 1500 -s mafft-nexus-internal-trimmed-gblocks1-clean-50p-IQTree.phylip --prefix iqtree-GHOST-50p -m GTR+FO*H4 -T AUTO --rcluster 10 --mrate G,R,E
```
For the older version (v. 1), use this command:

```
iqtree -st DNA -ninit 10 -bb 1500 -s mafft-nexus-edge-trimmed-gblocks-clean-50p-IQTree.phylip -pre iqtree-GHOST-50p -m GTR+FO*H4 -rcluster 10 -mrate G,R,E
```

### 4.2 ExaBayes
To infer the tree using Bayesian statistics, we will use the ExaBayes program. This program requires a significant amount of memory and is a slow process. Therefore, a compiler is used to perform parallel processes (MPI).

Additionally, we will need to download the [config.nex](https://github.com/uriborrajo/HETGEN1000/blob/main/config.nex) file, which specifies the parameters and variables for our Bayesian analysis. Review this file before using it and adjust it according to your analysis requirements.

To automate the process, we have created a [script](https://github.com/uriborrajo/HETGEN1000/blob/main/exabayes.sh) that performs four executions successively. We also recommend using the ```tmux``` program or another program to manage background processes, as this is a time-consuming process (if you review the [exabayes.sh](https://github.com/uriborrajo/HETGEN1000/blob/main/exabayes.sh) script, you will see the elapsed execution time for our analysis).

To run the script, execute the following command:

```
./exabayes.sh alignments concatenated_matrix.phylip config.nex
```

When run has completed, calculate these in an interactive session:
 
- Check if parameters converged (ESS should be >100, PSRF should be <1.1)
``` 
postProcParam -f ExaBayes_parameters.17taxa.0 ExaBayes_parameters.17taxa.1 -n 17params
```
- Calculate standard deviation of split frequencies (aka convergence, ASDSF should be <1%)
```
sdsf -f ExaBayes_topologies.run-0.24Ony-6 ExaBayes_topologies.run-1.24Ony-6
```
- Compute credible set of topologies (aka frequency of diff topologies)
```
credibleSet -f ExaBayes_topologies.17taxa.0 ExaBayes_topologies.17taxa.1 -n 17Cred
```
- Extract bipartitions (see how well supported different nodes are)
```
extractBips -f ExaBayes_topologies.17taxa.0 ExaBayes_topologies.17taxa.1 -n 17Bips
```
*Output files to care about: .bipartitions = name of node; .bipartitionStatistics = support for node (check ESS)*
- Generate consensus tree
```
consense –f ExaBayes_topologies.16Taxa.0 ExaBayes_topologies.16Taxa.1 –n 16ConsTree
```

### 4.3 ASTRAL
```
iqtree2 -S mafft-nexus-internal-trimmed-gblocks1-clean-50p/ --prefix loci -T AUTO --seqtype DNA -m GTR+FO*H4 -B 1500
```

```
java -jar /home/intern/Desktop/apps/ASTRAL/astral.5.7.8.jar -i *.treefile -o astral_sptree.treefile
```

## 5. REFERENCES

• Crotty, S.M., Minh, B.Q., Bean, N.G., Holland, B.R., Tuke, J., Jermiin, L.S., Haeseler, A.V., 2019. GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol. 69, 249–264. https://doi.org/10.1093/sysbio/syz051

• Wang, H.C., Minh, B.Q., Susko, S., Roger, A.J., 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. Syst. Biol., 67:216-235. https://doi.org/10.1093/sysbio/syx068

• Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., Haeseler, A.V., Jermiin, L.S., 2017. ModelFinder: Fast Model Selection for Accurate Phyloge- netic Estimates, Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285

• Faircloth, B.C., 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. https://doi.org/10.1093/bioinformatics/btv646.

• Andre J. Aberer, Kassian Kobert, Alexandros Stamatakis, ExaBayes: Massively Parallel Bayesian Tree Inference for the Whole-Genome Era, Molecular Biology and Evolution, Volume 31, Issue 10, October 2014, Pages 2553–2556, https://doi.org/10.1093/molbev/msu236

• Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, Bui Quang Minh, IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies, Molecular Biology and Evolution, Volume 32, Issue 1, January 2015, Pages 268–274, https://doi.org/10.1093/molbev/msu300

• Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics (Oxford, England), 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565










































