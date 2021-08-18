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
	parser.add_argument('--min', '-n', required=False, help='Min read length - for MultiSeg')
	parser.add_argument('--total', '-t', required=True, help='Total reads')
	parser.add_argument('--copybackratio', '-cbr', required=False, help='copyback ratio; 5 copyback, 5 snapback, 3 copyback, 3 snapback (*comma separated, must add up to 1)  - for Copyback')

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
	
	if args.min:
		minLength = int(args.min)
	else:
		minLength = 300 #Default 300 nt minimum

	if args.total:
		totalReads = int(args.total)
	else:
		print("Total number of reads")

	if args.copybackratio:
		ratioList = [args.copybackratio]
		if sum(ratioList) == 1:
			cb5 = ratioList[0]
			sb5 = ratioList[1]
			cb3 = ratioList[2]
			# No need to have sb3 as all add up to 1
		else:
			print("Copy back ratio MUST add up to 1; 5 copyback, 5 snapback, 3 copyback, 3 snapback")
	else:
		#Default parameters
		cb5 = 0.45
		sb5 = 0.05
		cb3 = 0.45
		

	# Send to methods
	if args.simMethod:
		method = args.simMethod
		if method == "ViReMa":
			Simulate.ViReMa(refGenome, outputDir, maxLength, totalReads)
		if method == "INDEL":
			Simulate.INDEL(refGenome, outputDir, maxLength, totalReads)
		if method == "CopyBack":
			Simulate.CopyBack(refGenome, outputDir, maxLength, totalReads, cb5, sb5, cb3)
		if method == "MultiSeg":
			Simulate.MultiSeg(refGenome, outputDir, maxLength, minLength, totalReads)
	else:
		print("Must have method, either:\nViReMa\nINDEL\nCopyBack\nMultiSeg")

if __name__ == "__main__":
	main(sys.argv[1:])