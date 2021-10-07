#!/usr/bin/env python3

import argparse, sys

def main(argv):
	parser = argparse.ArgumentParser(description="Program to take results from DI-Tecor and create bed file")

	parser.add_argument("-d", required=True, help="DI-Tector Output Directory")
	
	args = parser.parse_args()

	if args.d:
		makeBed(args.d)
	else:
		print("Requires DI-Tector Output Directory")

def makeBed(dir):
	ditFile = open((dir + "/DI-tector_counts.txt"), "r")

	outF = open((dir + "/DI-tector.bed"), "w")
	outF.write("track name=Virus_Recombinations description=""Virus_Recombinations"" graphType=junctions\n")

	for line in ditFile:
		if "=" in line or line == "" or "Length" in line: # Skipping lines without DVG info
			continue

		lineList = line.split("\t")

		segList = lineList[5]
		segList = segList.split("|")
		if segList[0] != segList[1]: # Must be on the same segment
			continue
		if "Rev" in lineList[0]:
			strand = "-"
		else:
			strand = "+"

		outF.write(segList[0] + "\t" + lineList[2] + "\t" + lineList[3] + "\tNAMES_TBD\t" + lineList[6] + "\t" + strand + "\t" + lineList[2] + "\t" + lineList[3] + "\n")

	outF.close()
	print("Bed file saved in DI-Tector directory")


if __name__ == "__main__":
	main(sys.argv[1:]) 
