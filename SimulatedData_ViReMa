import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ViReMa simulation code 


fastaFile = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

n = int(maxLength/2)

reads = []
summary = []


for refGenome in SeqIO.parse(fastaFile, "fasta"):
	lenGenome = len(refGenome) - n # Get length of genome (# of nucleotides)
	print(lenGenome)

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


		print(type(dI))
		rec = SeqRecord(
    		Seq(dI,),
    		id="test|test|gb|ABC123.4|ABC123_4",
    		description="test DIPs",
		)	

		reads.append(rec)
		summary.append([count, pos1, pos1+n, pos2, pos2+n])

	SeqIO.write(reads, "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimulatedViReMa.fasta", "fasta")
	
	summaryFile = open("/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimulatedViReMa.csv", "w")

	count = 0
	while count < totalReads:
		print(count)
		summaryFile.write(str(summary[count][0]) + "\t" + str(summary[count][1]) + "\t" + str(summary[count][2]) + "\t" + str(summary[count][3]) + "\t" + str(summary[count][4]) + "\n")
		count += 1
	summaryFile.close()
