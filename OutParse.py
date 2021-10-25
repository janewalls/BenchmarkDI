#!/usr/bin/env python3

import argparse, sys


 

def main(argv):
	
	parser = argparse.ArgumentParser(description="A program to generate Defective Interfering particles with 4 different methods")

	# Create flag arguments
	parser.add_argument("-d", "--ditector", required=True, help="DI-Tector output file: DI-tector_output_sorted.txt")
	parser.add_argument("-v", "--virema", required=True, help="ViReMa output file: Virus_Recombination_Results.txt")
	parser.add_argument("-o", "--outdir", required=True, help="Output Directory")
	parser.add_argument("-s", "--sim", required=False, help="Simualated data")
	parser.add_argument("-f", "--fasta", required=False, help="Simulated reads fasta file")
	parser.add_argument("--keep", required=False, help="Keep DIPs where tools say BP and RI are the same in ViReMa, default off")
	parser.add_argument("--std", required=False, type=int, help="Standard Deviation")
	parser.add_argument("--cut", required=False, type=int, help="Cut off if want to remove DIPs that have low occurances, default = 0")	

	args = parser.parse_args()

	dataDict = {} # Initiate nested dictionary for info

	totalList = []
	
	# Assign flag inputs to variables
	if args.keep:
		keep = True
	else:
		keep = False

	if args.ditector:
		outDIT = args.ditector
		ditParse(dataDict, outDIT, totalList) # Puts DI-Tector results in dictionary
	if args.virema:
		outVir = args.virema
		virParse(dataDict, outVir, totalList, keep) # Puts ViReMa results in dictionary
	if args.std:
		std = args.std
	else:
		std = 0
	if args.cut:
		cutOff = args.cut
	else:
		cutOff = 0
	if args.fasta:
		fasta = args.fasta
	else:
		fasta = "none"
	if args.outdir:
		parsed_out = args.outdir
	else:
		print("Requires output directory!")
	if args.sim:
		actual = args.sim
		simParse(dataDict, actual, totalList) # Puts simulated DIP info in dictionary
		compile(parsed_out, dataDict, std, fasta, totalList) # Compares tools with simulated data
	else:
		compare(parsed_out, dataDict, std, cutOff, totalList) # Compares tools between themselves
	




def simParse(dataDict, simData, totalList): # Takes csv and puts into searchable dictionary
	simDict = {}
	simRef= {} # Contains read numbers for later fasta
	simData = open(simData, "r")
	simData.readline() # Skips header line

	totCount = 0 # Total DIP read counter

	for line in simData: # reading each line
		lineList = line.split("\t") # creating a list of items in each line
		totCount += 1

		# Finds break point / ri and their respective segments
			# Differing locations depending on simulation type
		if len(lineList) == 7:
			bp = int(lineList[3])
			ri = int(lineList[5])
			segA = int(lineList[1])
			segB = int(lineList[4])
		elif len(lineList) == 4:
			bp = int(lineList[2])
			ri = int(lineList[3])
			segA = int(lineList[1])
			segB = int(lineList[1])
		else:
			bp = int(lineList[3])
			ri = int(lineList[4])
			segA = int(lineList[1])
			segB =	int(lineList[1])		
	
		# Puts information as key in dictionary
		if (segA,bp,segB,ri) in simDict:
			simDict[(segA,bp,segB,ri)] +=1 # If already exists adds count
		else:
			simDict[(segA,bp,segB,ri)] = 1 # If hasn't been seen
			simRef[(segA,bp,segB,ri)] = lineList[0]

	# Update dictionarys and lists
	dataDict.update({"simDict": simDict}) 
	dataDict.update({"simRef": simRef})
	totalList.append(["simDict", totCount])



def ditParse(dataDict, outDIT, totalList): # Parses DI-tector output
	ditDict = {}
	try:
		openDIT = open(outDIT, "r")
	except FileNotFoundError:
		print("File not found or DI-Tctor found no DIPs: check DI-Tector run")
		dataDict.update({"ditDict": ditDict})
		totalList.append(["ditDict", 0])
		return
	openDIT.readline()

	totCount = 0

	for line in openDIT:
		lineList = line.split("\t") # Splits line into list
		bp = int(lineList[2])
		ri = int(lineList[3])
		segA = int(lineList[8].split("seg")[1])
		segB = int(lineList[9].split("seg")[1])

		totCount += 1
		# Puts information on to dictionary with info as key
		if (segA, bp, segB, ri) in ditDict:
			ditDict[(segA, bp, segB, ri)] += 1 # If already seen
		else:
			ditDict[(segA, bp, segB, ri)] = 1 # If not seen, creates new item

	dataDict.update({"ditDict": ditDict})
	totalList.append(["ditDict", totCount])


def virParse(dataDict, outVir, totalList, keep): # Parses ViReMa output
	virDict = {}
	outVir = open(outVir, "r")

	totCount = 0

	for line in outVir:
		# Find out which segment
		if line[:1] == "@": 
			if "@EndofLibrary" in line:
				continue
			lineList = line.split(" ")
			lineList = lineList[1].split("_")
			# Exctract segment information
			segA = int(lineList[0].split("seg")[1])
			if lineList[1] == "RevStrand":
				#segA += "R"
				segB = int(lineList[3].split("seg")[1])
				#if len(lineList) == 4:
					#segB += "R"
					
			else:
				segB = int(lineList[2].split("seg")[1])
				#if len(lineList) == 3:
					#segB += "R"

		#Extract DIPs	
		else:
			lineList = line.split("\t")
			lineList.pop()
			for i in lineList: # Final item in list is a "/n"
				item = i.split("_")
				bp = int(item[0])
				ri = int(item[2])
				if keep == False and segA == segB and bp == ri:
					continue # Ensures DIPs found in virema that have identical bp and ri positions will be removed
				virDict[(segA, bp, segB, ri)] = int(item[4])
				totCount += int(item[4])

	dataDict.update({"virDict": virDict})
	totalList.append(["virDict", totCount])


def compile(outputPath, dataDict, sd, fasta, totalList): # puts into CSV alongside simulated reads

	# Open file and write header 
	output = open((outputPath + "/parser_output.csv"), "w")
	output.write("BP_seg" + "\t" + "BP" + "\t" + "RI_seg" + "\t" + "RI" + "\t" + "Sim_Count" + "\t" + "DITector_count" + "\t" + "ViReMa_count" + "\n")
	
	simDict = dataDict["simDict"]
	ditDict = dataDict["ditDict"]	
	virDict = dataDict["virDict"]

	unidList = []

	# Total read, and DIP matches
	ditTotMatch = 0
	virTotMatch = 0
	ditDipMatch = 0
	virDipMatch = 0

	ditUnct = 0
	ditDipUnct = 0
	
	for key in  simDict.keys(): # Checks all simulated DIPs
		ditCount = 0
		virCount = 0

		# Checks if ViReMa found DIP, within STD - searches dictionary
		filtDict= {k: v for k, v in virDict.items() if (k[0] == key[0] and k[2] == key[2]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1]) and (k[3]+sd >= key[3] and k[3]-sd <= key[3])}
		for x in filtDict:
			virCount += filtDict[x]
			virTotMatch += filtDict[x]
			virDipMatch +=1
			del virDict[x]  # delete item in list to ensure not double counted
		
		# Checks if ViReMa found DIP, within STD - searches dictionary
		filtDict= {k: v for k, v in ditDict.items() if (k[0] == key[0] and k[2] == key[2]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1]) and (k[3]+sd >= key[3] and k[3]-sd <= key[3])}
		for x in filtDict:
			ditCount += filtDict[x]
			ditTotMatch += filtDict[x]
			ditDipMatch += 1
			del ditDict[x]  # Same as above

		# Write DIP line in parser output file
		outString = str(key[0]) + "\t" + str(key[1]) + "\t" + str(key[2]) + "\t" + str(key[3]) + "\t" + str(simDict[key]) + "\t" + str(ditCount) + "\t" + str(virCount) + "\n"
		output.write(outString)

		# If not found in either tool
		if ditCount == 0 and virCount == 0:
			unidList.append([key])

	output.close()

	# Creates fasta file of unfound DIPs
	if fasta != "none":
		fasta = open(fasta, "r")
		outFasta = open((output + "/unfound_fasta.csv"), "w")
		simRef = dataDict["simRef"]
		reads = fasta.split(">")
		for key in unidList:
			pos = simRef[key]

			outFasta.write(">" + reads[pos] + "\n")
	outFasta.close()

	# Creates summary file 
	summaryFile = open((outputPath + "/Summary_results.txt"), "w")

	# Finds ViReMa Unmatched DIPs/Reads
	virUnct = 0	
	virDipUnct = 0	
	for x in virDict:
		virUnct += virDict[x]
		virDipUnct += 1

	# Finds DI-Tector Unmatched DIPs/Reads
	ditUnct = 0
	ditDipUnct = 0
	for x in ditDict:
		ditUnct += ditDict[x]
		ditDipUnct += 1
	
	# Creating strings for summary file 
	totSim = "Total sim reads" + "\t" + str(totalList[2][1]) + "\n"
	totDit = "Total DITector Reads" + "\t" + str(totalList[0][1]) + "\n"
	totVir = "Total ViReMa Reads" + "\t" + str(totalList[1][1]) + "\n"

	matDit = "DITector Matched Reads" + "\t" + str(ditTotMatch) + "\n"
	matVir = "ViReMa Matched Reads" + "\t" + str(virTotMatch) + "\n"
	unDit = "DITector UnMatched Reads" + "\t" + str(ditUnct) + "\n"
	unVir = "ViReMa UnMatched Reads" + "\t" + str(virUnct) + "\n"

	simDip = "Total DIPs" + "\t" + str(len(simDict.keys())) + "\n"
	matDipDit = "DITector Matched DIPs" + "\t" + str(ditDipMatch) + "\n"
	matDipVir = "ViReMa Matched DIPs" + "\t" + str(virDipMatch) + "\n"
	unDipDit = "DITector UnMatched DIPs" + "\t" + str(ditDipUnct) + "\n"
	unDipVir = "ViReMa UnMatched DIPs" + "\t" + str(virDipUnct) 
	
	# Write summary file
	summaryFile.write(totSim + totDit + totVir + matDit + matVir + unDit + unVir + simDip + matDipDit + matDipVir + unDipDit + unDipVir)
	summaryFile.close()

	print("\nAction Complete, File Saved\n")


def compare(outputPath, dataDict, sd, cutOff, totalList): # puts into CSV without simulated data

	# Open file and write header
	output = open((outputPath + "/parser_output.csv"), "w")
	output.write("BP_seg" + "\t" + "BP" + "\t" + "RI_seg" + "\t" + "RI" + "\t" + "DITector_count" + "\t" + "ViReMa_count" + "\n")
	
	ditDict = dataDict["ditDict"]	
	virDict = dataDict["virDict"]

	# Number of DIPs in each DIP from the respective tools
	totDipDit = len(ditDict) 
	totDipVir = len(virDict)

	totMatch = 0
	dipMatch = 0

	delList = []

	for key in ditDict.keys(): # Runs though DI-Tector DIPs
		virCount = 0
		# Checks if ViReMa found as well - searches dictionary
		filtDict= {k: v for k, v in virDict.items() if (k[0] == key[0] and k[2] == key[2]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1]) and (k[3]+sd >= key[3] and k[3]-sd <= key[3])}
		for x in filtDict:
			virCount += filtDict[x] 
			totMatch += filtDict[x]
			dipMatch +=1
			del virDict[x] # delete item in list to ensure not double counted
		if virCount > cutOff or ditDict[key] > cutOff: # Notes DIP if count/occurance above cut off
			outString = str(key[0]) + "\t" + str(key[1]) + "\t" + str(key[2]) + "\t" + str(key[3]) + "\t" + str(ditDict[key]) + "\t" + str(virCount) + "\n"
			output.write(outString)
		delList.append(key)
	
	# Deletes items that have been matched with ViReMa
	for key in delList:
		del ditDict[key]

	
	virUnct = 0	
	virDipUnct = 0	
	for x in virDict:
		virUnct += virDict[x]
		virDipUnct += 1
		# Adds extra ViReMa counts that aren't found in DI-Tector, so long as cut off allows
		if cutOff == 0:
			outString = str(x[0]) + "\t" + str(x[1]) + "\t" + str(x[2]) + "\t" + str(x[3]) + "\t" + str(0) + "\t" + str(virDict[x]) + "\n"
			output.write(outString)

	output.close()

	# Find number of DIPs that didn't match with ViReMa
	ditUnct = 0
	ditDipUnct = 0
	for x in ditDict:
		ditUnct += ditDict[x]
		ditDipUnct += 1

	summaryFile = open((outputPath + "/Summary_results.txt"), "w")

	# Creating strings for summary file 
	totDit = "Total DITector Reads" + "\t" + str(totalList[0][1]) + "\n"
	totVir = "Total ViReMa Reads" + "\t" + str(totalList[1][1]) + "\n"
	totDipDit = "Total DITector DIPs" + "\t" + str(totDipDit) + "\n"
	totDipVir = "Total ViReMa DIPs" + "\t" + str(totDipVir) + "\n"

	mat = "Matched Reads" + "\t" + str(totMatch) + "\n"
	unDit = "UnMatched Dip Reads" + "\t" + str(ditUnct) + "\n"
	unVir = "UnMatched Vir Reads" + "\t" + str(virUnct) + "\n"

	matDip = "Matched DIPs" + "\t" + str(dipMatch) + "\n"
	unDipDit = "DITector UnMatched DIPs" + "\t" + str(ditDipUnct) + "\n"
	unDipVir = "ViReMa UnMatched DIPs" + "\t" + str(virDipUnct)

	# Write summary file
	summaryFile.write(totDit + totVir + totDipDit + totDipVir + mat + unDit + unVir + matDip + unDipDit + unDipVir)
	summaryFile.close()

	print("\nAction Complete, File Saved\n")
	

if __name__ == "__main__":
	main(sys.argv[1:]) 

