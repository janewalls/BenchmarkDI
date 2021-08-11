import Bio, random

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def ViReMa(refGenome, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)
    
	reads = []
	readOutput = outDir + "/SimViReMa.fasta"
	summaryOutput = outDir + "/SimViReMa.csv"

	summaryFile = open(summaryOutput, "w")

	n = int(maxLen/2)

	count = 0

	# Create reads

	lenGenome = len(refGenome) - n # Get length of genome (# of nucleotides)
	assNo = refGenome.id

	while count < totalReads:
		while True:
			pos1 = random.randint(0,lenGenome)
			pos2 = random.randint(0,lenGenome)

			if ((pos2 > pos1 + n) or (pos1 > pos2 + n)):
				break

		seq1 = str(refGenome.seq[pos1:(pos1+n)])
		seq2 = str(refGenome.seq[pos2:(pos2+n)])

		count += 1
		dI = seq1 + seq2

        
		rec = SeqRecord(
			Seq(dI,),
			id=assNo,
			description=str(pos1) + "-" + str(pos1+n) + "," + str(pos2) + "-" + str(pos2+n),
		)	

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(pos1) + "\t" + str(pos1+n) + "\t" + str(pos2) + "\t" + str(pos2+n) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 
    
	summaryFile.close()


def INDEL(refGenome, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)

 
	reads = []
	readOutput = outDir + "/SimDITectorIndel.fasta"
	summaryOutput = outDir + "/SimDITectorIndel.csv"

	summaryFile = open(summaryOutput, "w")

	count = 0

	#Create Reads
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id

	while count < totalReads:
		while True:
			frag1 = random.randint(1,maxLen)
			frag2 = maxLen - frag1
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
				id=assNo,
				description=str(bP-frag1) + "-" + str(bP) + "," + str(rI) + "-" + str(rI+frag2),
			)	

			reads.append(rec) # Adds new reads to list

			summaryFile.write(str(count) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

	summaryFile.close()


def CopyBack(refGenome, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)


	reads = []
	readOutput = outDir + "/SimDITectorCopyBack.fasta"
	summaryOutput = outDir + "/SimDITectorCopyBack.csv"

	summaryFile = open(summaryOutput, "w")

	count = 0
    
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id

    # 5' Copy Backs (45% of total reads)
	while count < int(totalReads*0.45):
		while True:
			j = random.randint(1,maxLen)
			frag1 = int(j/2)
			frag2 = maxLen - j
			bP = random.randint(frag1,(lenGenome - frag2))
			rI = random.randint(frag1,(lenGenome - frag2))
			
			if bP < rI:
				break
			
		seq1 = refGenome.seq[(bP-frag1):bP] # Finds sequence
		seq2 = str(seq1.reverse_complement()) # reverse compliment on 5' end
		seq3 = str(refGenome.seq[rI:(rI+frag2)])

		count += 1
		dI = str(seq1) + seq3 + seq2

		rec = SeqRecord(
			Seq(dI,),
			id=assNo,
			description=str(bP-frag1) + "-" + str(bP) + "," + str(rI) + "-" + str(rI+frag2) + "," + str(bP) + "-" + str(bP-frag1),
		)

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\t" + str(bP) + "\t" + str(bP-frag1) + "\n") # Save summary info to file

    # 5' Snap back (5% of total reads)
	while count < int(totalReads*0.5):
		j = random.randint(1,maxLen)
		frag = random.randint(1,(lenGenome - j))
			
		seq1 = refGenome.seq[frag:(frag+j)]
		seq2 = str(seq1.reverse_complement()) # reverse compliment on 5' end

		count += 1
		dI = str(seq1) + seq2

		rec = SeqRecord(
			Seq(dI,),
			id=assNo,
			description=str(frag) + "-" + str(frag+j) + "," + str(frag+j) + "-" + str(frag),
		)

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(frag) + "\t" + str(frag+j) + "\t" + str(frag+j) + "\t" + str(frag) + "\n") # Save summary info to file

    # 3' Copy back (45% of total reads)
	while count < int(totalReads*0.95):
		while True:
			j = random.randint(1,maxLen)
			frag1 = int(j/2)
			frag2 = maxLen - j
			bP = random.randint(frag1,(lenGenome - frag2))
			rI = random.randint(frag1,(lenGenome - frag2))
			
			if bP < rI:
				break
			
		seq1 = str(refGenome.seq[(bP-frag1):bP]) 
		seq2 = refGenome.seq[rI:(rI+frag2)] # Finds sequence
		seq3 = str(seq2.reverse_complement()) # Then reverse compliment for 3'


		count += 1
		dI = seq3 + seq1 + str(seq2)

		rec = SeqRecord(
			Seq(dI,),
			id=assNo,
			description=str(rI) + "-" + str(rI+frag2) + "," + str(bP-frag1) + "-" + str(bP) + "," + str(rI+frag2) + "-" + str(rI),
		)

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI+frag2) + "\t" + str(rI) + "\n") # Save summary info to file


    # 3' Snap back (5% of total reads)
	while count < int(totalReads):
		j = random.randint(1,maxLen)
		frag = random.randint(1,(lenGenome - j))
			
		seq1 = refGenome.seq[frag:(frag+j)]
		seq2 = str(seq1.reverse_complement()) # reverse compliment on 5' end

		count += 1
		dI = seq2 + str(seq1)

		rec = SeqRecord(
			Seq(dI,),
			id=assNo,
			description=str(frag+j) + "-" + str(frag) + "," + str(frag) + "-" + str(frag+j),
		)

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(frag+j) + "\t" + str(frag) + "\t" + str(frag) + "\t" + str(frag+j) + "\n") # Save summary info to file


	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

	summaryFile.close()



def MultiSeg(records, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)

	reads = []
	readOutput = outDir + "/SimMultiSeg.fasta"
	summaryOutput = outDir + "/SimMultiSeg.csv"

	summaryFile = open(summaryOutput, "w")	

	count = 0

	while count < totalReads:
		ranSeg1 = random.randint(0,len(records)-1)
		seg1 = records[ranSeg1]

		bP = random.randint(0,600)

		ranSeg2 = random.randint(0,len(records)-1)
		seg2 = records[ranSeg2]

		rI = random.randint(len(seg2) - 600,len(seg2))


		seq1 = str(seg1.seq[1:bP])
		seq2 = str(seg2.seq[rI:len(seg2)])

		count += 1
		dI = seq1 + seq2

		rec = SeqRecord(
			Seq(dI,),
			id=str(records[ranSeg1]),
			description="seg" + str(ranSeg1+1) + ":1-" + str(bP) + "," + str(ranSeg2+1) + str(rI) + "-" + str(len(seg2)),
		)	

		reads.append(rec) # Adds new reads to list
		summaryFile.write(str(count) + "\t" + str(ranSeg1+1) + "\t" + str(1) + "\t" + str(bP) + "\t" + str(ranSeg2+1) + "\t" + str(rI) + "\t" + str(len(seg2)) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 
        
	summaryFile.close()