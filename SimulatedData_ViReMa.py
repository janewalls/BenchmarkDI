import random, sys, getopt
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ViReMa simulation code 

def main(argv):

	# Initialise adjustable parameters 
	fastaFile = ""
	outputDir = ""
	maxLength = 0
	totalReads = 0
	count = 0



	# Get options/paramaters 
	try:
		opts, args = getopt.getopt(argv,"h:i:o:m:t:",["ifile=","odir=","mlen=", "tread="])
	except getopt.GetoptError:
		print("script.py -i <inputfile> -o <outputdir> -m <maxreadlength> -t <totalreads>")
		sys.exit(1)
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

	# Initialise variables necessary for files
	reads = []
	readOutput = outputDir + "/SimulatedViReMa.fasta"
	summaryOutput = outputDir + "/SimulatedViReMa.csv"

	summaryFile = open(summaryOutput, "w")

	n = int(maxLength/2)

	# Create reads
	for refGenome in SeqIO.parse(fastaFile, "fasta"):
		lenGenome = len(refGenome) - n # Get length of genome (# of nucleotides)

		while count <= totalReads:
			while True:
				pos1 = random.randint(1,lenGenome)
				pos2 = random.randint(1,lenGenome)

				if ((pos2 > pos1 + n) or (pos1 > pos2 + n)):
					break

			seq1 = str(refGenome.seq[pos1:(pos1+n)])
			seq2 = str(refGenome.seq[pos2:(pos2+n)])

			count += 1
			dI = seq1 + seq2


			rec = SeqRecord(
				Seq(dI,),
				id="test|test|gb|ABC123.4|ABC123_4",
				description="test DIPs",
			)	

			reads.append(rec)

			summaryFile.write(str(count) + "\t" + str(pos1) + "\t" + str(pos1+n) + "\t" + str(pos2) + "\t" + str(pos2+n) + "\n")

		SeqIO.write(reads, readOutput, "fasta")
		
		summaryFile.close()

if __name__ == "__main__":
	main(sys.argv[1:])