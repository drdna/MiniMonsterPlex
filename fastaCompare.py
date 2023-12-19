sampleData = {}
with open("AllDataReorderedAnnotatedNR.fasta", 'r') as read:
	gate = False
	for line in read:
		if gate:
			newLine = line.replace('\n', '')
			sampleData[oldLine.split('_')[0]] = newLine
			gate = False

		if line[0:3] == '>UF':
			oldLine = line
			gate = True

myData = []
with open("builtSeq.fasta", 'r') as read:
	gate = False
	for line in read:
		if gate:
			myData.append([oldLine, line.replace('\n', '')])
			gate = False

		if line[0:3] == '>UF':
			oldLine = line.replace('_', "").split('h')[0]
			gate = True

#print(sampleData)

#print(myData)
with open("fastaCompare.txt", 'a') as append:
	for sample in myData:
		if sample[0] in sampleData:
			append.write(sample[0] + '\n')
			append.write("my data:" + sample[1] + '\n')
			append.write("sample :" + sampleData[sample[0]] + '\n')
