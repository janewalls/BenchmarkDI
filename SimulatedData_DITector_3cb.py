import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# DITector - 5' copy back simulation code 

# Adjustable parameters 
fastaFile = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

readOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.fasta"
summaryOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.csv"

#Initialising lists
reads = []

summaryFile = open(summaryOutput, "w")

for refGenome in SeqIO.parse(fastaFile, "fasta"):
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)

	while count <= totalReads:
		while True:
			frag1 = random.randint(1,maxLength)
			frag2 = maxLength - frag1
			bP = random.randint(frag1,(lenGenome - frag2))
			rI = random.randint(frag1,(lenGenome - frag2))
			
			if bP < rI:
				break
		
		seq1 = str(refGenome.seq[(bP-frag1):bP]) 

		seq2 = str(refGenome.seq[rI:(rI+frag2)].reverse_compliment()) # Finds sequence, then reverse compliment

		count += 1
		dI = seq1 + seq2

		rec = SeqRecord(
			Seq(dI,),
			id="test|test|gb|ABC123.4|ABC123_4",
			description="test DIPs",
		)	

		reads.append(rec)

		summaryFile.write(str(count) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\n")

	SeqIO.write(reads, readOutput, "fasta")

	summaryFile.close()
