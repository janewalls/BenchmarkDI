import argparse, sys, Simulate
from os import error
from Bio import SeqIO

def main(argv):

	parser = argparse.ArgumentParser()

	# Positional 
	parser.add_argument('simMethod', help='Method for simulation')

	# Flagged arguments
	parser.add_argument('--file', '-f', required=True, help='Input fasta')
	parser.add_argument('--outdir', '-o', required=True, help='Output Directory')
	parser.add_argument('--max', '-m', required=True, help='Max read length')
	parser.add_argument('--total', '-t', required=True, help='Total reads')

	args = parser.parse_args()

	# Parsing arguments and assigning to variables
	if args.file:
		fastaFile = args.file
		try:
			refGenome = SeqIO.read(fastaFile, "fasta")
		except ValueError:
			refGenome = list(SeqIO.parse(fastaFile, "fasta"))
	else:
		print("Fasta file needed!")

	if args.outdir:
		outputDir = args.outdir
	else:
		print("Output directory needed")

	if args.max:
		maxLength = args.max
	else:
		print("Max read length")

	if args.total:
		totalReads = args.total
	else:
		print("Total number of reads")

	# Send to methods
	if args.simMethod:
		method = args.simMethod
		if method == "ViReMa":
			Simulate.ViReMa(refGenome, outputDir, maxLength, totalReads)
		if method == "INDEL":
			Simulate.INDEL(refGenome, outputDir, maxLength, totalReads)
		if method == "CopyBack":
			Simulate.CopyBack(refGenome, outputDir, maxLength, totalReads)
		if method == "MultiSeg":
			Simulate.MultiSeg(refGenome, outputDir, maxLength, totalReads)
	else:
		print("Must have method, either:\nViReMa\nINDEL\nCopyBack\nMultiSeg")

if __name__ == "__main__":
	main(sys.argv[1:])