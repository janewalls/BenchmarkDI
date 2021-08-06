import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ViReMa simulation code 

# Adjustable parameters 
fastaFile = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

readOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimulatedViReMa.fasta"
summaryOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimulatedViReMa.csv"

# Initialise 
reads = []

summaryFile = open(summaryOutput, "w")

n = int(maxLength/2)


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
