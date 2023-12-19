import os
import glob

fileList = glob.glob('fastq/*.gz')

for file in fileList:
	fileNum = file.split('/')[1].split('.')[0]
	print(fileNum)
	os.system(f'hisat2 -p 2 -x 70-15index -U {file} --max-intronlen 20 --summary-file OutNew/{fileNum}summary.txt --dta-cufflinks | samtools sort - -@ 2 -O bam -o OutNew/{fileNum}hits.bam')
	os.system(f'tabix OutNew/{fileNum}hits.bam')
	os.system(f'bcftools mpileup --threads 2 -d 100000 -R MonsterPlexRegionsFileSuperCont.txt --annotate FORMAT/AD -f 70-15.fasta.fasta OutNew/{fileNum}hits.bam >> OutNew/{fileNum}.vcf')
	os.system(f'bcftools call -c --ploidy 1 OutNew/{fileNum}.vcf -o OutNew/seperateCall/{fileNum}call.vcf')
	os.system(f'bedtools genomecov -ibam OutNew/{fileNum}hits.bam -bg > OutNew/Coverage/{fileNum}cover.bed')
	os.system(f'bgzip OutNew/{fileNum}.vcf')
	os.system(f'tabix OutNew/{fileNum}.vcf.gz')
