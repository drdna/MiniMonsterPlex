import argparse

#setting up arg parsing for the output folder
parser = argparse.ArgumentParser(
    prog='FastaCompare',
    description=(
        'A tool for comparing fasta files'
    ),
    epilog=(
        'Arguments are required'
    ),
)

#command line option for input
parser.add_argument(
    '-i',
    action='store',
    nargs=2,
    help=(
        'The two fasta files to be compared.' 
        'This argument exepcts two file names seperated by a space.'
    ),
    required = True,
    default='n/a'
)

#command line option for output folder
parser.add_argument(
    '-o',
    action='store',
    help=(
        'User defined name for the output file.'
    ),
    default='fastaCompare.txt'
)

args = parser.parse_args()

output_File = args.o
file1 = args.i[0]
file2 = args.i[1]

sampleData = {}
with open(file1, 'r') as read:
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
with open(file2, 'r') as read:
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
with open(output_File, 'a') as append:
	for sample in myData:
		if sample[0] in sampleData:
			append.write(sample[0] + '\n')
			append.write("file1:" + sample[1] + '\n')
			append.write("file2:" + sampleData[sample[0]] + '\n')
