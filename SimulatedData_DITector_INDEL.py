import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# DITector - INDEL simulation code 

fastaFile = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

reads = []
summary = []


for refGenome in SeqIO.parse(fastaFile, "fasta"):
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	print(lenGenome)

	while count <= totalReads:
		while True:
			bP = random.randint(1,(lenGenome - maxLength))
            rI = random.randint(1,lenGenome)
			j = random.randint(j,maxLength)

			if ((ri > bp + j) or (bp > ri + j)):
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
		summary.append([count, bP, bP+j, (rI - (maxLength-j)), rI])

	SeqIO.write(reads, "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.fasta", "fasta")
	
	summaryFile = open("/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.csv", "w")

	count = 0
	while count < totalReads:
		print(count)
		summaryFile.write(str(summary[count][0]) + "\t" + str(summary[count][1]) + "\t" + str(summary[count][2]) + "\t" + str(summary[count][3]) + "\t" + str(summary[count][4]) + "\n")
		count += 1
	summaryFile.close()
