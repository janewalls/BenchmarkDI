import random, sys, getopt
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# DITector - INDEL simulation code 

def main(argv):

	# Initialise adjustable parameters 
	fastaFile = ""
	maxLength = 0
	totalReads = 0
	count = 0

# Get options/paramaters 
	try:
		opts, args = getopt.getopt(argv,"h:i:o:m:t:",["ifile=","odir=","mlen=", "tread="])
	except getopt.GetoptError:
		print("script.py -i <inputfile> -o <outputdir> -m <maxreadlength> -t <totalreads>")
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print("script.py -i <inputfile> -o <outputdir> -m <maxreadlength> -t <totalreads>")
			sys.exit()
		elif opt in ("-i", "--ifile"):
			fastaFile = str(arg)
		elif opt in ("-o", "--odir"):
			outputDir = str(arg)
		elif opt in ("-m", "--mlen"):
			maxLength = int(arg)
		elif opt in ("-t", "--tread"):
			totalReads = int(arg)


	#Initialising other
	reads = []
	readOutput = outputDir + "/SimulatedDITectorIndel.fasta"
	summaryOutput = outputDir + "/SimulatedDITectorIndel.csv"

	summaryFile = open(summaryOutput, "w")

	#Create Reads
	for refGenome in SeqIO.parse(fastaFile, "fasta"):
		lenGenome = len(refGenome) # Get length of genome (# of nucleotides)

		while count < totalReads:
			while True:
				frag1 = random.randint(1,maxLength)
				frag2 = maxLength - frag1
				bP = random.randint(frag1,(lenGenome - frag2))
				rI = random.randint(frag1,(lenGenome - frag2))
				
				if bP < rI:
					break
			
			seq1 = str(refGenome.seq[(bP-frag1):bP])
			seq2 = str(refGenome.seq[rI:(rI+frag2)])

			count += 1
			dI = seq1 + seq2

			rec = SeqRecord(
				Seq(dI,),
				id="test|test|gb|ABC123.4|ABC123_4",
				description="test DIPs",
			)	

			reads.append(rec) # Adds new reads to list

			summaryFile.write(str(count) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\n") # Save summary info to file

		SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

		summaryFile.close()

if __name__ == "__main__":
	main(sys.argv[1:])
