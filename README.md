# MiniMonsterPlex
MiniMonsterplex is an automatic variant calling pipeline. 

## Requirements 
Install via Conda:
* Python 3.6 or higher
* [Hisat2](https://anaconda.org/bioconda/hisat2)
* [Tabix](https://anaconda.org/bioconda/tabix)
* [Samtools](https://anaconda.org/bioconda/samtools)
* [Bcftools](https://anaconda.org/bioconda/bcftools)
* [BedTools](https://anaconda.org/bioconda/bedtools)

## Command Line Functionality 
```
Python3 MiniMonsterPlex.py -o [output folder name] -m [.csv metadata file name] -h
```
```-h```= Help command: including this flag will bring up the help screen.
```-o```= Output Folder: User given name for the created output folder. When no option is used it defaults to output. **Note** only two folders can have the same name so only one default folder can exist at a time.
```-m```= Metadata file: Name of the .csv metadata file formatted as shown below.

## Quick Start guide

Note Requires Python 3.6 or higher

Sample result files are contained within the folder sampleResults. These 4 files are all generated in this pipeline
1. Use the .gz files within the fastq folder.
   or (Optional) Extract the tarball WheatBlast.tar.gz into the fastq folder.
2. Run autoVCF.py: Produces coverage, bam, and vcf files.
   ```
   python3 autoVCF.py
   ```
3. Run autoMerge.py: produces tabix of all of the vcf files, a merged vcf file, and fastqListCall.txt a list of files fed into bcftools merge.
   ```
   python3 autoMerge.py
   ```
4. Run sampleBuilder.py: produces builtSeq.fasta and MonsterPlexSitesUsed.txt (a testing file used to compare the sites anayalzed compared to those within the original files).
   SampleBuilder now has support for metadata addition to fasta headers. This is via the .tsv file 'metadata.tsv'. Dummy metadata is included as an example. Make sure the sequenceID given in the metadata is the exact same as the file name without the extension.
   ```
   python3 sampleBuilder.py
   ```
6. (Optional) Run fastaCompare.py if you wish to see the built seq compared with sample results: produces fastaCompare.txt which overlays both sequences in an easy to parse format.
   ```
   python3 fastaCompare.py
   ```
