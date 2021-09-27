#!/usr/bin/env python3

import argparse, sys


#/Users/janewalls/Documents/VS_CODE/MastersProjectDI/BenchmarkDI/BenchmarkDI/OutParse.py -s "/Users/janewalls/OneDrive - University of Glasgow/Project/testParser/test1.csv" -d "/Users/janewalls/Documents/DI-tector_output_sorted.txt" -v "/Users/janewalls/Documents/Virus_Recombination_Results.txt" -o "/Users/janewalls/Documents" --std 2



def main(argv):
	
	parser = argparse.ArgumentParser(description="A program to generate Defective Interfering particles with 4 different methods")

	parser.add_argument("-s", "--sim", required=False, help="Simualated data")
	parser.add_argument("-d", "--ditector", required=True, help="DI-Tector output file: DI-tector_output_sorted.txt")
	parser.add_argument("-v", "--virema", required=True, help="ViReMa output file: Virus_Recombination_Results.txt")
	parser.add_argument("-o", "--outdir", required=True, help="Output Directory")
	parser.add_argument("-f", "--fasta", required=False, help="Simulated reads fasta file")
	parser.add_argument("--std", required=False, type=int, help="Standard Deviation")
	parser.add_argument("--cut", required=False, type=int, help="Cut off if want to remove DIPs that have low occurances")	

	args = parser.parse_args()

	dataDict = {} # Initiate nested dictionary

	totalList = []

	if args.ditector:
		outDIT = args.ditector
		ditParse(dataDict, outDIT, totalList)
	if args.virema:
		outVir = args.virema
		virParse(dataDict, outVir, totalList)
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
		simParse(dataDict, actual, totalList)

		compile(parsed_out, dataDict, std, fasta, totalList)
	else:
		compare(parsed_out, dataDict, std, cutOff, totalList)
	



# bp ri simCount, ditCount, viCount

def simParse(dataDict, simData, totalList): # puts into CSV alongside simulated reads
	simDict = {}
	simRef= {} # Contains read numbers for later fasta
	simData = open(simData, "r")
	simData.readline()

	totCount = 0 # Total DIP read count

	for line in simData:
		lineList = line.split("\t")
		totCount += 1

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
	
		if (segA,bp,segB,ri) in simDict:
			simDict[(segA,bp,segB,ri)] +=1 
		else:
			simDict[(segA,bp,segB,ri)] = 1 
			simRef[(segA,bp,segB,ri)] = lineList[0]

	
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
		if (bp,ri) in ditDict:
			ditDict[(segA, bp, segB, ri)] +=1 # Will need to be different for Copyback, and fragment -- can either have first line of csv to show method, adjust so all are the same, or put flag
		else:
			ditDict[(segA, bp, segB, ri)] = 1 

	dataDict.update({"ditDict": ditDict})
	totalList.append(["ditDict", totCount])


def virParse(dataDict, outVir, totalList): # Parses DI-tector output
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
					#print()
		#Extract DIPs	
		else:
			lineList = line.split("\t")
			lineList.pop()
			for i in lineList: # Final item in list is a "/n"
				item = i.split("_")
				bp = int(item[0])
				ri = int(item[2])
				virDict[(segA, bp, segB, ri)] = int(item[4])
				totCount += int(item[4])

	dataDict.update({"virDict": virDict})
	totalList.append(["virDict", totCount])


def compile(outputPath, dataDict, sd, fasta, totalList): # puts into CSV alongside simulated reads

	output = open((outputPath + "/parser_output.csv"), "w")
	output.write("Simulated Read (seg#,bp,seg#,ri)" + "\t" + "Simulated count" + "\t" + "DI-Tector count" + "\t" + "ViReMa count" + "\n" )
	
	simDict = dataDict["simDict"]
	ditDict = dataDict["ditDict"]	
	virDict = dataDict["virDict"]

	unidList = []

	ditTotMatch = 0
	virTotMatch = 0
	ditDipMatch = 0
	virDipMatch = 0

	ditUnct = 0
	ditDipUnct = 0
	
	for key in  simDict.keys():
		ditCount = 0
		virCount = 0

		filtDict= {k: v for k, v in virDict.items() if (k[0]+sd >= key[0] and k[0]-sd <= key[0]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1])}
		for x in filtDict:
			virCount += filtDict[x]
			virTotMatch += filtDict[x]
			virDipMatch +=1
			del virDict[x] # Same as above
		
		filtDict= {k: v for k, v in ditDict.items() if (k[0]+sd >= key[0] and k[0]-sd <= key[0]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1])}
		for x in filtDict:
			ditCount += filtDict[x]
			ditTotMatch += filtDict[x]
			ditDipMatch += 1
			del ditDict[x]  # items will only be counted once if sd > 0, and will be counted under first shown


		outString = str(key) + "\t" + str(simDict[key]) + "\t" + str(ditCount) + "\t" + str(virCount) + "\n"
		output.write(outString)

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

	# Creates summary file 
	summaryFile = open((outputPath + "/Summary_results.txt"), "w")

	virUnct = 0	
	virDipUnct = 0	
	for x in virDict:
		virUnct += virDict[x]
		virDipUnct += 1

	ditUnct = 0
	ditDipUnct = 0
	for x in ditDict:
		ditUnct += ditDict[x]
		ditDipUnct += 1
	
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
	
	summaryFile.write(totSim + totDit + totVir + matDit + matVir + unDit + unVir + simDip + matDipDit + matDipVir + unDipDit + unDipVir)
	summaryFile.close()

	

	print("\nAction Complete, File Saved\n")


def compare(outputPath, dataDict, sd, cutOff, totalList): # puts into CSV alongside simulated reads

	output = open((outputPath + "/parser_output.csv"), "w")
	output.write("Read (seg#,bp,seg#,ri)" + "\t" + "DI-Tector count" + "\t" + "ViReMa count" + "\n")
	
	ditDict = dataDict["ditDict"]	
	virDict = dataDict["virDict"]

	totDipDit = len(ditDict)
	totDipVir = len(virDict)

	totMatch = 0
	dipMatch = 0

	delDict = {}

	for key in  ditDict.keys():
		virCount = 0
		filtDict= {k: v for k, v in virDict.items() if (k[0]+sd >= key[0] and k[0]-sd <= key[0]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1])}
		for x in filtDict:
			virCount += filtDict[x]
			totMatch += filtDict[x]
			dipMatch +=1
			del virDict[x] # Same as above
		if virCount > cutOff or ditDict[key] > cutOff:
			outString = str(key) + "\t" + str(ditDict[key]) + "\t" + str(ditDict[key]) + "\t" + str(virCount) + "\n"
			output.write(outString)
		delDict[key] = ditDict[key]
	
	for key in delDict:
		del ditDict[key]

	virUnct = 0	
	virDipUnct = 0	
	for x in virDict:
		virUnct += virDict[x]
		virDipUnct += 1

	ditUnct = 0
	ditDipUnct = 0
	for x in ditDict:
		ditUnct += ditDict[x]
		ditDipUnct += 1

	summaryFile = open((outputPath + "/Summary_results.txt"), "w")

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

	summaryFile.write(totDit + totVir + totDipDit + totDipVir + totDipVir + mat + unDit + unVir + matDip + unDipDit + unDipVir)
	summaryFile.close()

if __name__ == "__main__":
	main(sys.argv[1:]) 

