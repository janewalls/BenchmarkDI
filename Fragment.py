import random 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

fastaFile="/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimMultiSeg.fasta"

fragOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/fragmentTest.fasta"

allReads = list(SeqIO.parse(fastaFile, "fasta"))

n = 10 # number of fragments per read

fragments = []

for read in allReads:
    for i in range(n):
        l = random.randint(300,400)
        if l > len(read): 
            break
        pos = random.randint(0,(len(read)-l)) #Finds start position for read
        frag = str(read.seq[pos:(pos+l)]) # Gets sequence

        rec = SeqRecord(
			Seq(frag,),
			id=str(read.id),
			description=str(pos) + "-" + str(pos+l),
		)	

        fragments.append(rec)

SeqIO.write(fragments, fragOutput, "fasta") # Save fragments to fasta file 


