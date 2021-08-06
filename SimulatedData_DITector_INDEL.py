import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# DITector - INDEL simulation code 

# Adjustable parameters 
fastaFile = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

readOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.fasta"
summaryOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.csv"

#Initialising lists
reads = []
summary = []

summaryFile = open(summaryOutput, "w")

for refGenome in SeqIO.parse(fastaFile, "fasta"):
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	print(lenGenome)

	while count <= totalReads:
		while True:
			j = random.randint(1,maxLength)			
			bP = random.randint(1,(lenGenome - maxLength))
			rI = random.randint(maxLength,lenGenome)
			
			if ((rI - (maxLength-j) > bP + j) or (bP > rI)):
				break
		
		seq1 = str(refGenome.seq[bP:(bP+j)])
		seq2 = str(refGenome.seq[(rI - (maxLength-j)):rI])

		count += 1
		dI = seq1 + seq2


		print(type(dI))
		rec = SeqRecord(
			Seq(dI,),
			id="test|test|gb|ABC123.4|ABC123_4",
			description="test DIPs",
		)	

		reads.append(rec)

		summaryFile.write(count + "\t" + bP + "\t" + bP+j + "\t" + (rI - (maxLength-j)) + "\t" + rI + "\n")

	SeqIO.write(reads, readOutput, "fasta")

	summaryFile.close()
