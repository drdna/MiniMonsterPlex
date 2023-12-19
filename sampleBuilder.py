import os

sites =[]
sitesUsed =[]
with open('MonsterPlexSitesList.txt', 'r') as read:
	for line in read:
		sites.append(line.strip('\n'))

#for site in sites:
#	if site == 'supercont8.5 2841082':
#		print('site 2841082 found')
#		quit()

with open('OutNew/seperateCall/wheatBlastMergedCallAll.vcf', 'r') as read:
	seqs = list()
	check = False
	for line in read:
		if line.split('\t')[0] == '#CHROM':
			print("header past seqs made")
			fqList = line.split('\t')
			for n in range(9,len(fqList)):
				seqs.append([fqList[n], ''])
			check = True
		elif check and line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1] in sites:
			#this creates a horizontal split of the line
			lineList = line.strip('\n').split('\t')
			for n in range(9,len(lineList)):
				fields = lineList[n].split(':')
#				if line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1] == 'supercont8.5 2841082':
#					print(fields)
#					quit()
				#this checks if no reads were found the locus
				if fields[0] == '.':
					seqs[n - 9][1] += "N"
				#this checks if reads for only ref or alt were found
				elif len(fields[2].split(',')) == 1:
					#this checks ref
					if fields[0] == '0':
						if int(fields[2]) > 5:
							seqs[n - 9][1] += lineList[3]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					#this checks alt
					elif fields[0] == '1':
						if int(fields[2]) > 5:
							seqs[n - 9][1] += lineList[4]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					else:
						print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
						print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						quit()
				#this checks cases were both ref and alt are registered
				elif len(fields[2].split(',')) >= 2:
					#this creates a list out of the AD field
					AD = fields[2].split(',')
					#this checks if ref is blank
					if AD[0] == '.':
						if int(AD[1]) > 5:
							seqs[n - 9][1] += lineList[4].split(',')[0]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					#this checks if alt is blank
					elif AD[1] == '.':
						if int(AD[0]) > 5:
							seqs[n - 9][1] += lineList[3]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					#checks if ref is greater then alt
					elif int(AD[0]) > int(AD[1]):
						if int(AD[0]) > (int(AD[1]) * 20):
							seqs[n - 9][1] += lineList[3]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					#checks if alt is greater than ref
					elif int(AD[1]) > int(AD[0]):
						if int(AD[1]) > (int(AD[0]) * 20):
							seqs[n - 9][1] += lineList[4].split(',')[0]
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						else:
							seqs[n - 9][1] += "N"
							sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					elif int(AD[1]) == int(AD[0]):
						seqs[n - 9][1] += lineList[3]
						sitesUsed.append(line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					else:
						print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
						print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
						quit()
				else:
					print("something went wrong with " + seqs[n][0] + '\n' + lineList[n])
					print("at site " + line.strip('\n').split('\t')[0] + ' ' + line.strip('\n').split('\t')[1])
					quit()

metaData = {}
with open('metadata.tsv', 'r') as read:
	for line in read:
		ID = line.split('\t')[0].strip('\n')
		species = line.split('\t')[1].strip('\n')
		host = line.split('\t')[2].strip('\n')
		country = line.split('\t')[3].strip('\n')
		metaData[ID] = [species, host, country]
		print(metaData[ID][0])

with open('MonsterPlexSitesUsed.txt', 'a') as write:
	old = 'l'
	for x in sitesUsed:
		if x != old:
			write.write(x + '\n')
			old = x
			print(x)

with open('builtSeqMeta.fasta', 'a') as writeSeq:
	for read in seqs:
		seqID = read[0].split('/')[1].split('.')[0]
		if seqID in metaData:
			seqSpecies = metaData[seqID][0]
			seqHost = metaData[seqID][1]
			seqCountry = metaData[seqID][2]
			writeSeq.write(f'>{seqID}_{seqSpecies}_{seqHost}_{seqCountry}\n{read[1]}\n')
		else:

			writeSeq.write('>' + read[0].split('/')[1].split('.')[0] + '\n' + read[1] + '\n')
print(seqs)
