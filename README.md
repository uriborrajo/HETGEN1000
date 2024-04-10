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

Phyluce is a software package initially designed to analyse data extracted from ultra-conserved elements in the genomes of organisms. These ultraconserved elements (UCE) are highly conserved regions of the genomes of organisms, in our case organisms within the phylum Mollusca.

For the present study we will use this software package to analyse data from 35 mollusc species.

To use phyluce you need to have MinicondaX installed (in our case Miniconda 3):
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```
If you want to install another version of miniconda, just change the name of the .sh in the wget command (e.g. Miniconda3 to Miniconda2).

Next, you must initialise Miniconda3:
```
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
To verify that Miniconda is installed and working properly, close and reopen your terminal and then run ```conda list```, this should produce an output.

The installation of Miniconda3 is only necessary to be able to separate the downloaded packages in different environments, avoiding that they are all mixed in the base environment. Thus, Conda allows us to create different environments where we will install Phyluce in an environment called ```Phyluce-1.7.2``` where ```-1.7.2``` is the installed version.

To install phyluce-1.7.2 and create the new environment:

```
wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.2/distrib/phyluce-1.7.2-py36-Linux-conda.yml
conda env create -n phyluce-1.7.2 --file phyluce-1.7.2-py36-Linux-conda.yml
```

Once the software package is installed, we are ready to activate and use the corresponding environment. It should be noted that whenever we want to use phyluce software we will have to activate the environment where it is installed.

To activate the environment:

```
conda activate phyluce-1.7.2
```
or
```
source activate phyluce-1.7.2
```

### 2. IQ-TREE INSTALLATION
An inference of phylogenetic trees by the maximum likelihood method is one of the programs that IQ-TREE is often used for. The program offers an accurate and time efficient scholastic algorithm, the same accuracy and time effectiveness are shared with other phylogenetic inference programs, e.g., such as RAxML and PhyML (Nguyen et al., 2015).

To install IQ-TREE we will use the conda command, as it allows us to download this software program with a simple command:
```
conda install -c bioconda iqtree
```
### 3. EXABAYES INSTALLATION
ExaBayes is a bioinformatics program for Bayesian phylogenetic analysis. It was especially inspired by MrBayes, but also applies similar approaches from BEAST. It sets up a Markov chain Monte Carlo sampling method that allows to determine parameters of the evolutionary model (e.g. branch length or substitution rates) and the posterior probability of a tree or topology (Andre J. et al., 2014).

For the installation of ExaBayes follow the installation instructions (section 3) on the official ExaBayes website by clicking on this link: https://cme.h-its.org/exelixis/web/software/exabayes/manual/manual.html#sec-11

It should be noted that the installation of ExaBayes did not work with versions higher than gcc version 10 and we had to downgrade it. 

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
CD-HIT is a program that applies an algorithm to reduce size redundancy and thereby enhance the performance of advanced sequence analysis methods. It is based on the prediction of a similarity to filter out the extraneous conversion of sequences to their debugging counterparts (Fu L. et al., 2012).

For CD-HIT installation:
```
conda install -c bioconda cd-hit
```
For CD-HIT-DUP and CD-HIT-EST installation: 
```
conda install -c bioconda cd-hit-auxtools
```
For the present study, only CD-HIT-DUP has been used. This is used before making assemblies with Spades.

>CD-HIT is used after making assemblies with Spades.

>CD-HIT-EST is used for protein, i.e. transcriptome alignments.


### 5. ASTRAL INSTALLATION
The ASTRAL tool is a java-based package that is used for estimating a phylogenetic tree with a set of unrooted gene-trees. It detect the species tree that turns out to be the one which gives the majority of induced quartet trees within present gene trees, and this is a statistical approach which is compatible with a multispecies coalescent model (visit Astral's GitHub because a new version called ASTER has been released).

To install Astral, we create a folder called Apps, which is going to be our central place for all the extensions and programs we need during the course.

Then: 

```
mkdir -p Apps/Astral
cd Apps/Astral
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip
```
It should be noted that whenever we want to use Astral, we will have to specify the path where the Astral directory with the executable script of the program is located (later you will see in the commands that the path to the Astral directory is used).

### 6. FASTP INSTALLATION

```
conda install -c bioconda fastp
```


### 3. COUNT THE READ DATA
- To count the number of reads in a sequence file for a species, Unix tools can be used. 

- This command counts reads from R1, which should be equal to R2. Otherwise, we will have to equal the total reads of the two sequences (see UNPAIRED FASTQ FILES). 

- The command divides the output number by 4 to get the number of sequence reads.

```
for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
```

### 4. FASTP

To use fastp we will have to run the script ```fastp.sh```, this script uses as input the folder with the fastqs. Also, the output can be called clean-fastq to find out where the cleaned sequences are.

An example for running the script would be:
```
./fastp.sh ~/Desktop/Oriol/fastq ~/Desktop/Oriol/clean-fastq
```
Note that if you copy the script directly into a new shell script, you will have to give it permissions so that it can be executed. To give it these permissions:
```
chmod +x fastp.sh
```

### 5. CD-HIT-DUP
Before making the Spades assemblies we have to clean possible contaminations and duplicates of the raw data.

· Run the cd-hit-dup.sh shell script as
```
./cd-hit-dup.sh raw-data cdhitdup
```

###### FASTQ COMBINE PAIRED END (UNPAIRED FASTQ FILES)
For those unpaired fastq files (i.e. those files where R1 and R2 do not have the same number of reads) we will have to compare and discard the read-only data.

IMPORTANT - You have to verify that the fastq files use exactly four lines per sequence, otherwise this program will not recognise these sequences.

The program's output will be three files, the first two contains the reads of the matching pairs and the third one the reads of the matching pairs.

Usage: python fastqCombinePairedEnd.py input1 input2 separator

input1 = LEFT  fastq or fastq.gz file (R1)

input2 = RIGHT fastq or fastq.gz file (R2)

separator

```
while read i; do \
cd /home/intern/Desktop/data/UCE_clean_reads/${i}/split-adapter-quality-trimmed; \
python ../../../scripts/Oriol/fastqCombinePairedEnd.py ${i}-READ1.fastq ${i}-READ2.fastq separator; \
done < especies_sin_archivos.txt
```

### 6. SPADES

Steps to crate "assembly.conf" file:
```
cd cdhitdup
```

``` 
echo "[samples]" > ../assembly.conf
```

``` 
for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/"; done >> ../assembly.conf
```

or

``` 
for i in *; do echo "$i:intern/home/Desktop/data/cdhitdup/$i/split-adapter-quality-trimmed/"; done >> ../assembly.conf
```


```
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --memory 20000 \
    --cores 30 \
```

>### WARNING!!
>
>Based on [Brant's](https://gist.github.com/brantfaircloth/e48e7e4eb9748854962863d104f94095) python script we kept UCE loci that match more than one contig.
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
>```
>mv duplicates.fasta ../taxon-sets/all
>```
>```
>phyluce_assembly_explode_get_fastas_file \
>    --input duplicates.fasta \
>    --output exploded-fastas \
>    --by-taxon
>```
>```
>cd exploded-fastas
>```
>```
>cat *-DUPE1.unaligned.fasta >> duplicates.fasta
>```
>```
>sed -i 's/_DUPE1//g' duplicates.fasta
>```
>```
>cd ../
>```
>```
>cat all-taxa-incomplete.fasta duplicates.fasta >> all-taxa-incomplete2.fasta
>```


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
We are going to trim these loci using **Gblocks**

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

![NEW](https://img.shields.io/badge/New-blue?style=social)
### ZORRO

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
















## REFERENCES

• Crotty, S.M., Minh, B.Q., Bean, N.G., Holland, B.R., Tuke, J., Jermiin, L.S., Haeseler, A.V., 2019. GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol. 69, 249–264. https://doi.org/10.1093/sysbio/syz051

• Wang, H.C., Minh, B.Q., Susko, S., Roger, A.J., 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. Syst. Biol., 67:216-235. https://doi.org/10.1093/sysbio/syx068

• Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., Haeseler, A.V., Jermiin, L.S., 2017. ModelFinder: Fast Model Selection for Accurate Phyloge- netic Estimates, Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285

• Faircloth, B.C., 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. https://doi.org/10.1093/bioinformatics/btv646.

• Andre J. Aberer, Kassian Kobert, Alexandros Stamatakis, ExaBayes: Massively Parallel Bayesian Tree Inference for the Whole-Genome Era, Molecular Biology and Evolution, Volume 31, Issue 10, October 2014, Pages 2553–2556, https://doi.org/10.1093/molbev/msu236

• Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, Bui Quang Minh, IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies, Molecular Biology and Evolution, Volume 32, Issue 1, January 2015, Pages 268–274, https://doi.org/10.1093/molbev/msu300

• Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics (Oxford, England), 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565










































