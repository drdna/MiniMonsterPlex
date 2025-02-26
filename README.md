# MiniMonsterPlex
MiniMonsterplex is an automatic variant calling pipeline. 

## Table of Contents
1. [Requirements](https://github.com/TrStans606/MiniMonsterPlex/blob/main/README.md#requirements)
2. [Data Input](https://github.com/TrStans606/MiniMonsterPlex/blob/main/README.md#data-input)
3. [Command Line Functions](https://github.com/TrStans606/MiniMonsterPlex/blob/main/README.md#command-line-functions)
4. [Metadata Format](https://github.com/TrStans606/MiniMonsterPlex/tree/main#metadata-format)
5. [RAXML Help](https://github.com/TrStans606/MiniMonsterPlex#raxml)
6. [Tree building with MLtree](https://github.com/TrStans606/MiniMonsterPlex/tree/main#treebuilding-with-mltree)

## Requirements 
Install via Conda:
* Python 3.6 or higher
* R 3.2.1 or higher
* R package: [ape](https://cran.r-project.org/web/packages/ape/index.html)
* R package: [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
* R Bioconductor package: [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
* R package: [ggtext](https://cran.r-project.org/web/packages/ggtext/index.html)
* [Hisat2](https://anaconda.org/bioconda/hisat2)
* [Tabix](https://anaconda.org/bioconda/tabix)
* [Samtools](https://anaconda.org/bioconda/samtools)
* [Bcftools](https://anaconda.org/bioconda/bcftools)
* [BedTools](https://anaconda.org/bioconda/bedtools)
Build from source:
* [Standard RAXML Sequential version:](https://github.com/stamatak/standard-RAxML)

## Data Input
Fastq files with either a .fq or .fastq extension should be gzip compressed, extension .gz, and dropped into the [fastq/](fastq) folder before running. If your files are all uncompressed try using this command in the [fastq/](fastq) folder to compress them in bulk:
```
bgzip *.fastq
```
or
```
bgzip *.fq
```
Depending on what extension your files are.

## Command Line Functions
```
python3 MiniMonsterPlex.py -o [output folder name] -m [.csv metadata file name] -r [raxml binary name] -i [isolate_1] [isolate_2] -il [example.txt] -hf [host_1] [host_2] -hfl [example.txt] -h
```
+ ```-h```= Help command: including this flag will bring up the help screen.
+ ```-o```= Output Folder: User given name for the created output folder. When no option is used it defaults to output. **Note** only two folders can have the same name so only one default folder can exist at a time.
+ ```-m```= Metadata file: Name of the .csv metadata file formatted as shown below.
+ ```-r```=Raxml version: the name of the standard raxml binary

Filtering options:

**RAxML requires a minimum of four isolates in a multi-fasta file to generate a tree. If you do not provide four isolates or your selected host population(s) does not have four isolates the program will stop and ask if you want to continue without filtering or quit entirely.**

NOTE: Isolate should be the name of the file you are uploading minus the extenesions: so SRR1571.fq.gz will be SRR1571. Host names should be the exact same as those entered into your metadata file.

+ ```-i```= Isolate list[Optional]: a space separated list of all isolates you want included in the tree building. 
+ ```-il```= Isolate file[Optional]: a newline separated .txt file of all isolates you want included in the tree building. This can be combined with -i.
+ ```-hf```= Host list[Optional]: a space separated list that includes the hosts for the isolates you want included in tree building.
+ ```-hfl```= Host file[Optional]: a newline separated .txt file that lists the hosts of the isolates you want included in the tree building. This can be combined with -hf.

The host and isolate filtering can be combined. In that case the program will first filter by host and then filter by isolate. 

## Metadata Format

MiniMonsterPlex requires a custom .csv format for metadata:
```
sampleID,species,host,country
104,Po,Oryza,China
```
* ```sampleID``` should be the exact prefix of the fastq file given to MiniMonsterPlex so in the provided example it would be *104.fastq*.
* ```species``` is the name of the species to which the sample belongs.
* ```host``` is the host genus of the specific isolate.
* ```country``` is the country of origin.

A sample csv file can be found here [metadata.csv](metadata.csv)

## RAxML 
* Build the Standard Sequential RAXML version listed in [Requirements](https://github.com/TrStans606/MiniMonsterPlex/blob/main/README.md#requirements) from source
* The Standard, SSE3, or AVX version will work.
* Move the created exectuable into the MiniMonsterPlex directory
* Provide the exact name of the RAxML binary as the -r argument

## TreeBuilding with MLtree

![mlTree_sample](https://github.com/TrStans606/MiniMonsterPlex/assets/100236022/f6d01b13-eb93-42f3-80e8-d21ade5a5689)
