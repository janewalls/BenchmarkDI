#!/usr/bin/env python3

import argparse, sys, Generate
import re
from Bio import SeqIO

def main(argv):

	parser = argparse.ArgumentParser(description="A program to generate Defective Interfering particles with 4 different methods")

	# Positional 
	parser.add_argument("simMethod", help="Method for simulation; either ViReMa, INDEL, Copyback, MultiSeg")

	# Flagged arguments
	parser.add_argument("-f", "--file", required=True, help="Input fasta")
	parser.add_argument("-o", "--outdir", required=True, help="Output Directory")
	parser.add_argument("-m", "--max", type=int, required=True, help="Max read length")
	parser.add_argument("-t", "--total", type=int, required=True, help="Total reads")
	
	parser.add_argument("-c", "--copybackratio", required=False, help="copyback ratio; 5 copyback, 5 snapback, 3 copyback, 3 snapback (*comma separated, must add up to 1)  - for Copyback")
	parser.add_argument("--seg", type=int, required=False, help="For single segment use: fasta file number, only use if multisegment fasta file is input")
	parser.add_argument("--min", "-n", type=int, required=False, help="Min read length - for MultiSeg")
	parser.add_argument("--fragment", action="store_true", required=False, help="MultiSeg only - Fragment reads")
	parser.add_argument("-x", "--num", type=int, required=False, help="MultiSeg only - Number of fragments, default = 100000")
	parser.add_argument("-l", "--len", type=int, required=False, help="MultiSeg only - Average read length, default = 300")
	parser.add_argument("-s", "--std", type=int, required=False, help="MultiSeg only - Length standard deviation, default = 50")

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
		minLength = args.min
	else:
		minLength = 300 #Default 300 nt minimum

	if args.total:
		totalReads = args.total
	else:
		print("Total number of reads")

	if args.copybackratio:
		ratioList = args.copybackratio.split(",")
		ratioList = [float(i) for i in ratioList] # Cast list to integer
		if sum(ratioList) == 1:
			cb5 = ratioList[0]
			sb5 = ratioList[1]
			cb3 = ratioList[2]
			# No need to have sb3 as all add up to 1
		else:
			print("Copy back ratio MUST add up to 1; 5 copyback,5 snapback,3 copyback,3 snapback")
	else:
		#Default parameters
		cb5 = 0.45
		sb5 = 0.05
		cb3 = 0.45

	if args.seg:
		seg = args.seg
		refGenome = refGenome[seg-1]
		print("DIPs produced from segment "+ refGenome.id)
	else:
		while type(refGenome) == list:
			if args.simMethod == "MultiSeg":
				break # Allows for multi segments for option
			refGenome = refGenome[0]

	if args.fragment:
		frag = True
	else:
		frag = False
	
	if args.num:
		num = args.num
	else:
		num = 100000
	
	if args.len:
		len = args.len
	else:
		len = 300

	if args.std:
		std = args.std
	else:
		std = 50

	# Send to methods
	if args.simMethod:
		method = args.simMethod
		if method == "ViReMa":
			Generate.ViReMa(refGenome, outputDir, maxLength, totalReads)
		if method == "INDEL":
			Generate.INDEL(refGenome, outputDir, maxLength, totalReads)
		if method == "CopyBack":
			Generate.CopyBack(refGenome, outputDir, maxLength, totalReads, cb5, sb5, cb3)
		if method == "MultiSeg":
			Generate.MultiSeg(refGenome, outputDir, maxLength, minLength, totalReads, frag, num, len, std)
	else:
		print("Must have method, either:\nViReMa\nINDEL\nCopyBack\nMultiSeg")

if __name__ == "__main__":
	main(sys.argv[1:])