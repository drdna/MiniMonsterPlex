# MiniMonsterPlex
Repository for files used in MiniMonsterPlex pipeline

## Quick Start guide

Note Requires Python 3.6 or higher

Sample result files are contained within the folder sampleResults. These 4 files are all generated in this pipeline

1. Extract the tarball WheatBlast.tar.gz into the fastq folder.
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
