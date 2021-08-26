import random, argparse
from typing_extensions import Required
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Adjustable parameters 
fastaFile = ""
outputDir = ""
totalReads = 0
count = 0

readOutput = ""
summaryOutput = ""

reads = []


parser = argparse.ArgumentParser()

parser.add_argument('--file', '-f', required=True, help='Input fasta')
parser.add_argument('--outdir', '-o', required=True, help='Output Directory')
parser.add_argument('-m', required=True, help='Max length')
parser.add_argument('-t', required=True, help='Total reads')

args = parser.parse_args()

if args.file:
	fastaFile = args.file
else:
    print("Fasta file needed!")

if args.outdir:
	outputDir = args.outDir
else:
    print("Output directory needed")

for refGenome in SeqIO.parse(fastaFile, "fasta"):
	records = list(refGenome)

	assNo = refGenome.id

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
			id=assNo,
			description=str("seg" + (ranSeg1+1) + ":1," + bP + "-" + rI + "," + len(seg2)),
		)	

		reads.append(rec) # Adds new reads to list

		summaryOutput.write(str(count) + "\t" + str(1) + "\t" + str(bP) + "\t" + str(rI) + "\t" + str(len(seg2)) + "\n") # Save summary info to file

	SeqIO.write(reads, readOutput, "fasta") # Save reads to fasta file 
    
	summaryOutput.close()
