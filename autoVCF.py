import subprocess
import glob
import os

fileList = glob.glob('fastq/*.gz')

for file in fileList:
    fileNum = file.split('/')[1].split('.')[0]
    print(fileNum)
    os.makedirs("Out/seperateCall/")
    os.makedirs("Out/Coverage/")
    os.system(f'hisat2 -p 2 -x 70-15index -U {file} --max-intronlen 20 --summary-file Out/{fileNum}summary.txt --dta-cufflinks | samtools sort - -@ 2 -O bam -o Out/{fileNum}hits.bam')
    os.system(f'tabix Out/{fileNum}hits.bam')
    os.system(f'bcftools mpileup --threads 2 -d 100000 -R MonsterPlexRegionsFileSuperCont.txt --annotate FORMAT/AD -f index/70-15.fasta.fasta Out/{fileNum}hits.bam >> Out/{fileNum}.vcf')
    os.system(f'bcftools call -c --ploidy 1 Out/{fileNum}.vcf -o Out/seperateCall/{fileNum}call.vcf')
    os.system(f'bedtools genomecov -ibam Out/{fileNum}hits.bam -bg > Out/Coverage/{fileNum}cover.bed')
    os.system(f'bgzip Out/{fileNum}.vcf')
    os.system(f'tabix Out/{fileNum}.vcf.gz')
