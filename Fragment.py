import random, sys, argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Initialise variables
fastaFile= ""
fragOutput = ""
n = 0 
l = 0
s = 0

fragments = []

def main(argv):

	parser = argparse.ArgumentParser()

	# Flagged arguments
	parser.add_argument('--file', '-f', required=True, help='Input fasta')
	parser.add_argument('--out', '-o', required=True, help='Output fasta')
	parser.add_argument('--num', '-n', required=True, help='Number of fragments per read')
	parser.add_argument('--len', '-l', required=False, help='Average read length')
	parser.add_argument('--std', '-s', required=True, help='Length standard deviation')

	args = parser.parse_args()

	# Parsing arguments and assigning to variables
	if args.file:
		fastaFile = args.file
	else:
		print("Fasta file needed!")

	if args.out:
		fragOutput = args.outdir
	else:
		print("Output fasta")

	if args.num:
		n = args.num
	else:
		print("Number of fragments per read needed")
	
	if args.len:
		l = args.len
	else:
		print("Average read length")

	if args.std:
		s = args.std
	else:
		print("Length standard deviation")


	allReads = list(SeqIO.parse(fastaFile, "fasta"))

	for read in allReads:
		for i in range(n):
			ln = random.randint(l-s,l+s)
			if ln > len(read): 
				break
			pos = random.randint(0,(len(read)-ln)) #Finds start position for read
			frag = str(read.seq[pos:(pos+ln)]) # Gets sequence

			rec = SeqRecord(
				Seq(frag,),
				id=str(read.id),
				description=str(pos) + "-" + str(pos+ln),
			)	

			fragments.append(rec)

	SeqIO.write(fragments, fragOutput, "fasta") # Save fragments to fasta file 

if __name__ == "__main__":
	main(sys.argv[1:])

