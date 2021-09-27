import random

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def MBP(refGenome, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)
    
	reads = []
	readOutput = outDir + "/SimMBP.fasta"
	summaryOutput = outDir + "/SimMBP.csv"

	summaryFile = open(summaryOutput, "w")

	summaryFile.write("Read #" + "\t" + "Seg#" + "\t" + "Start" + "\t" + "BP" + "\t" + "RI" + "\t" + "End" + "\n")

	n = int(maxLen/2)

	count = 0

	# Create reads

	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id
	seg = assNo.split("seg")[1]

	while count < totalReads:
		while True:
			pos1 = random.randint(0,(lenGenome - n))
			pos2 = random.randint(0,(lenGenome - n))

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

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(pos1) + "\t" + str(pos1+n) + "\t" + str(pos2) + "\t" + str(pos2+n) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 
    
	summaryFile.close()


def INDEL(refGenome, outDir, maxLen, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)
 
	reads = []
	readOutput = outDir + "/SimIndel.fasta"
	summaryOutput = outDir + "/SimIndel.csv"

	summaryFile = open(summaryOutput, "w")

	summaryFile.write("Read #" + "\t" + "Seg#" + "\t" + "Start" + "\t" + "BP" + "\t" + "RI" + "\t" + "End" + "\n")

	count = 0

	#Create Reads
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id
	seg = assNo.split("seg")[1]

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

			summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

	summaryFile.close()


def CopyBack(refGenome, outDir, maxLen, totalReads, cb5, sb5, cb3):
	maxLen = int(maxLen)
	totalReads = int(totalReads)


	reads = []
	readOutput = outDir + "/SimCopyBack.fasta"
	summaryOutput = outDir + "/SimCopyBack.csv"

	summaryFile = open(summaryOutput, "w")

	summaryFile.write("Read #" + "\t" + "Seg#" + "\t" + "Start" + "\t" + "BP" + "\t" + "RI" + "\t" + "loop_end" + "\t" + "BP" + "\t" + "End" + "\n")

	count = 0
    
	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id
	seg = assNo.split("seg")[1]

    # 5' Copy Backs (45% of total reads)
	while count < int(totalReads*cb5):
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

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\t" + str(bP) + "\t" + str(bP-frag1) + "\n") # Save summary info to file

    # 5' Snap back (5% of total reads)
	while count < int(totalReads*(sb5 +cb5)):
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

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(frag) + "\t" + str(frag+j) + "\t" + str(frag+j) + "\t" + str(frag) + "\t" + "-" + "\t" + "-" + "\n") # Save summary info to file

    # 3' Copy back (45% of total reads)
	while count < int(totalReads*(cb3 +sb5+cb5)):
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

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(rI+frag2) + "\t" + str(rI) + "\n") # Save summary info to file


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

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(frag+j) + "\t" + str(frag) + "\t" + str(frag) + "\t" + str(frag+j) + "\t" + "-" + "\t" + "-" + "\n") # Save summary info to file


	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

	summaryFile.close()



def MultiSeg(records, outDir, maxLen, minLen, totalReads, frag, fnum, flen, fstd):
	maxLen = int(maxLen)
	i = int(maxLen/2)
	totalReads = int(totalReads)

	reads = []
	readOutput = outDir + "/SimMultiSeg.fasta"
	summaryOutput = outDir + "/SimMultiSeg.csv"

	summaryFile = open(summaryOutput, "w")	

	summaryFile.write("Read#" + "\t" + "Seg#"  + "\t" + "Start" + "\t" + "BP" + "\t" + "Seg#" + "\t" + "RI" + "\t" + "End" + "\n")

	count = 0

	#for fragment
	bpList =[]

	while count < totalReads:
		while True:
			ranSeg1 = random.randint(0,len(records)-1)
			seg1 = records[ranSeg1]

			bP = random.randint(1,i)

			ranSeg2 = random.randint(0,len(records)-1)
			seg2 = records[ranSeg2]

			rI = random.randint(len(seg2) - i,len(seg2))

			if(bP + (len(seg2)-rI) >= minLen): # Ensures minimum threshold met
				break

		seq1 = seg1.seq[0:bP]
		seq2 = str(seg2.seq[rI:len(seg2)])

		count += 1
		dI = str(seq1) + seq2

		rec = SeqRecord(
			Seq(dI,),
			id=str(seg1.id + "," + seg2.id),
			description="seg" + str(ranSeg1+1) + ":1-" + str(bP) + "," + str(ranSeg2+1) + str(rI) + "-" + str(len(seg2)),
		)	

		bpList.append(bP)

		reads.append(rec) # Adds new reads to list
		summaryFile.write(str(count) + "\t" + str(ranSeg1+1) + "\t" + str(1) + "\t" + str(bP) + "\t" + str(ranSeg2+1) + "\t" + str(rI) + "\t" + str(len(seg2)) + "\n") # Save summary info to file


	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 
        
	summaryFile.close()

	if frag == True:
		Fragment(readOutput, bpList, outDir, fnum, flen, fstd)



def Fragment(fastaFile, bpList, outDir, num, ln, std):
	
	fragOutput = outDir + "/fragmented.fasta"
	summaryFile = outDir + "/fragmented.csv"

	summaryFile = open(summaryFile, "w")	

	records = list(SeqIO.parse(fastaFile, "fasta"))

	n = int(num/len(records)) # finds number of fragments per read needed

	fragments = []

	rCount = 0 # read counter
	fCount = 0 # fragment counter

	for read in records:
		bP = bpList[rCount]
		rCount += 1
		
		# creates n number of fragments per read
		for i in range(n): 
			rln = random.randint(ln-std,ln+std) # random fragment length
			while rln > len(read): 
				rln = random.randint(ln-std,ln+std) # In case fragment size is > read length

			pos = random.randint(0,(len(read)-rln)) #Finds start position for read
			frag = str(read.seq[pos:(pos+rln)]) # Gets sequence

			# Whether or not fragment spans over breakpoint
			if (bP > pos) & (bP < (pos+rln)):
				bPBool = True
			else:
				bPBool = False

			rec = SeqRecord(
				Seq(frag,),
				id=str(read.id),
				description=str(pos) + "-" + str(pos+rln),
			)	
			fCount += 1

			fragments.append(rec)
			summaryFile.write(str(rCount) + "\t" + str(fCount) + "\t" + str(rln) + "\t" + str(bPBool) + "\n") # Save summary info to file (read #, fragment #, fragment length, bp boolean)

	SeqIO.write(fragments, fragOutput, "fasta") # Save fragments to fasta file 

	summaryFile.close()


def MultiSeg2(records, outDir, maxLen, window, totalReads):
	maxLen = int(maxLen)
	totalReads = int(totalReads)

 
	reads = []
	readOutput = outDir + "/SimMultiSeg2.fasta"
	summaryOutput = outDir + "/SimMultiSeg2.csv"

	summaryFile = open(summaryOutput, "w")

	summaryFile.write("Read#" + "\t" + "Seg#"  + "\t" + "Start" + "\t" + "BP" + "\t" + "Seg#" + "\t" + "RI" + "\t" + "End" + "\n")

	count = 0


	while count < totalReads:
		while True:
			ranSeg1 = random.randint(0,len(records)-1)
			seg1 = records[ranSeg1]
			ranSeg2 = random.randint(0,len(records)-1)
			seg2 = records[ranSeg2]

			frag1 = random.randint(1,maxLen)
			frag2 = maxLen - frag1

			if window == 0:
				bP = random.randint(frag1,(len(seg1) - frag2))
				rI = random.randint(frag1,(len(seg2) - frag2))
			else:
				bP = random.randint(frag1,window)
				rI = random.randint((len(seg2)-window),(len(seg2)-frag2))

			if (ranSeg1 == ranSeg2 and bP > rI):
				break
			
			seq1 = seg1.seq[(bP-frag1):bP]
			seq2 = seg2.seq[rI:(rI+frag2)]

			count += 1
			dI = str(seq1) + str(seq2)

			rec = SeqRecord(
				Seq(dI,),
				id=str(seg1.id + "," + seg2.id),
				description=str(ranSeg1) + str(bP-frag1) + "-" + str(bP) + "," + str(ranSeg2) + str(rI) + "-" + str(rI+frag2),
			)	

			reads.append(rec) # Adds new reads to list

			summaryFile.write(str(count) + "\t" + str(ranSeg1) + "\t" + str(bP-frag1) + "\t" + str(bP) + "\t" + str(ranSeg2) + "\t" + str(rI) + "\t" + str(rI+frag2) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 

	summaryFile.close()


def NoDIP(refGenome, outDir, maxLen, totalReads):

	maxLen = int(maxLen)
	totalReads = int(totalReads)

	lenGenome = len(refGenome) # Get length of genome (# of nucleotides)
	assNo = refGenome.id
	seg = assNo.split("seg")[1]

	reads = []
	readOutput = outDir + "/SimNoDIP.fasta"
	summaryOutput = outDir + "/SimNoDIP.csv"

	summaryFile = open(summaryOutput, "w")
	summaryFile.write("Read#" "\t" + "Seg#" + "\t" + "Start" + "\t" + "End" + "\n")

	count = 0

	while count < totalReads:
		count += 1
		start = random.randint(0,(lenGenome - maxLen))
		seq = str(refGenome.seq[start:(start+maxLen)])
		
		rec = SeqRecord(
			Seq(seq,),
			id=assNo,
			description=str(start) + "-" + str(start+maxLen),
		)

		reads.append(rec) # Adds new reads to list

		summaryFile.write(str(count) + "\t" + str(seg) + "\t" + str(start) + "\t" + str(start+maxLen) + "\n") # Save summary info to file
	
	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file
	
	summaryFile.close()

