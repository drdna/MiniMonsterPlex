import os
import glob


fileList = glob.glob('fastq/*.gz')

for file in fileList:
	fileNum = file.split('/')[1].split('.')[0]
	print(fileNum)
	os.system(f'bgzip OutNew/seperateCall/{fileNum}call.vcf')
	os.system(f'tabix OutNew/seperateCall/{fileNum}call.vcf.gz')


with open('fastqListCall.txt', 'a') as append:
	for file in fileList:
		append.write(file.split('/')[1].split('.')[0] + 'call.vcf.gz\n')

os.system('bcftools merge -l fastqListCall.txt -o OutNew/seperateCall/wheatBlastMergedCallAll.vcf')
