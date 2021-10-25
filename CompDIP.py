#!/usr/bin/env python3

import argparse, sys
from collections import defaultdict

def main(argv):
	parser = argparse.ArgumentParser(description="A program to compare and compile difference sample outputs of Outparse.py: (5 samples max)")
	
	# Create flags 
	parser.add_argument("--s1", required=True, help="Sample parser_output.csv 1")
	parser.add_argument("--s2", required=True, help="Sample parser_output.csv 2")
	parser.add_argument("--s3", required=False, help="Sample parser_output.csv 3")
	parser.add_argument("--s4", required=False, help="Sample parser_output.csv 4")
	parser.add_argument("--s5", required=False, help="Sample parser_output.csv 5")
	parser.add_argument("-o", required=True, help="OutputFile")
	parser.add_argument("--cut", required=False, type=int, help="Cut off if want to remove DIPs that have low occurances, default to accept any occurance between either tool, default => 0")	

	args = parser.parse_args()

	sampleList = [] # Put samples in list, allows for variable number of samples

	# Assign flag user input to variables
	if args.s1:
		sampleList.append(args.s1)
	else:
		print("Requires an argument for at least s1 and s2")

	if args.s2:
		sampleList.append(args.s2)
	else:
		print("Requires an argument for at least s1 and s2")

	if args.s3:
		sampleList.append(args.s3)
	if args.s4:
		sampleList.append(args.s4)
	if args.s5:
		sampleList.append(args.s5)

	if args.o:
		outFile = args.o
	else:
		print("Requires output directory")

	if args.cut:
		cutOff = args.cut
	else:
		cutOff = 0
	compMethod(sampleList, outFile, cutOff) # Method that compiles samples

	print("Files compiled")


def compMethod(sList, outF, cut): 
	outF = open(outF, "w")
	header = "BP_seg\tBP\tRI_seg\tRI" # Create headers

	sDict = {}
	count=0

	for sample in sList: # Iterate through each sample file
		sample = open(sample,"r")

		sample.readline() # Read past header  

		sampleDict = {}
		sampleDict = defaultdict(list)

		for line in sample: # Reads individual samples
			lineList = line.split("\t")
			lineList = [int(i) for i in lineList]

			if (lineList[4]) >= cut or lineList[5] >= cut: # Filter out low occurances if set between tool
				sampleDict[(lineList[0],lineList[1],lineList[2],lineList[3])].append(lineList[4]) # DI-Tector Value, at DIP
				sampleDict[(lineList[0],lineList[1],lineList[2],lineList[3])].append(lineList[5]) # ViReMa Value, at DIP
				
		sDict.update({count: sampleDict}) # Add to all samples
		count += 1

		header += "\t" + "s" + str(count) + "_DIT" + "\t" + "s" + str(count) + "_VRM" # Adds appropriate file header for sample
	header += "\n"
	outF.write(header) # Write final heaer to file

	delList = [] #Inititate list

	# Compile all sample's dictionaries in file
	for sample in sDict: # iterate through sample
		
		# Deletes item in all samples in list from previous sample to ensure not repeated DIPs counted
		for i in sDict: 
			for item in delList:
				if item in sDict[i]:
					del sDict[i][item] # Ensure not to add dips twice 

		# Write DIP key
		for key in sDict[sample]:
			outF.write(str(key[0]) + "\t" + str(key[1]) + "\t" + str(key[2]) + "\t" + str(key[3]))
			
			# Iterate through sample and if DIP present notes count
			for i in sDict: 
				if key in sDict[i]:
					valList = sDict[i][key]
					outF.write("\t" + str(valList[0]) + "\t" + str(valList[1])) # Count for DI-Tector and ViReMa
				else:
					outF.write("\t" + str(0) + "\t" + str(0)) # If DIP not present, count will be 0 for both values
			
			delList.append(key) # Add key/DIP to delete list
			outF.write("\n")
	
	outF.close()

			


if __name__ == "__main__":
	main(sys.argv[1:]) 