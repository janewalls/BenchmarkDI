#!/usr/bin/env python3

import argparse, sys


#/Users/janewalls/Documents/VS_CODE/MastersProjectDI/BenchmarkDI/BenchmarkDI/OutParse.py -s "/Users/janewalls/OneDrive - University of Glasgow/Project/testParser/test1.csv" -d "/Users/janewalls/OneDrive - University of Glasgow/Project/sim2/SimViReMa_DIT/DI-tector_output_sorted.txt" -v "/Users/janewalls/OneDrive - University of Glasgow/Project/testParser/Virus_Recombination_Results.txt" -o "/Users/janewalls/OneDrive - University of Glasgow/Project/output.csv" --std 2


def main(argv):
	
	parser = argparse.ArgumentParser(description="A program to generate Defective Interfering particles with 4 different methods")

	parser.add_argument("-s", "--sim", required=False, help="Simualated data")
	parser.add_argument("-d", "--ditector", required=False, help="DI-Tector data")
	parser.add_argument("-v", "--virema", required=False, help="ViReMa")
	parser.add_argument("-o", "--outfile", required=True, help="Output Directory")
	parser.add_argument("--std", required=False, type=int, help="Standard Deviation")	

	args = parser.parse_args()

	dataDict = {} # Initiate nested dictionary

	if args.sim:
		actual = args.sim
		simParse(dataDict, actual)
	if args.ditector:
		outDIT = args.ditector
		ditParse(dataDict, outDIT)
	if args.virema:
		outVir = args.virema
		virParse(dataDict, outVir)
	if args.std:
		std = args.std
	if args.outfile:
		parsed_out = args.outfile
		compile(parsed_out, dataDict, std)
	



# bp ri simCount, ditCount, viCount

def simParse(dataDict, simData): # puts into CSV alongside simulated reads
	simDict = {}
	simData = open(simData, "r")

	for line in simData:
		lineList = line.split("\t")
		
		if len(lineList) == 7: # If using multiseg file
			bp = int(lineList[3])
			ri = int(lineList[5])
		else:
			bp = int(lineList[2])
			ri = int(lineList[3])
		
		if (bp,ri) in simDict:
			simDict[(bp,ri)] +=1 
		else:
			simDict[(bp,ri)] = 1 

	dataDict.update({"simDict": simDict})



def ditParse(dataDict, outDIT): # Parses DI-tector output
	ditDict = {}
	openDIT = open(outDIT, "r")
	openDIT.readline()

	for line in openDIT:
		lineList = line.split("\t") # Splits line into list
		bp = int(lineList[2])
		ri = int(lineList[3])

		if (bp,ri) in ditDict:
			ditDict[(bp,ri)] +=1 # Will need to be different for Copyback, and fragment -- can either have first line of csv to show method, adjust so all are the same, or put flag
		else:
			ditDict[(bp,ri)] = 1 

	dataDict.update({"ditDict": ditDict})


def virParse(dataDict, outVir): # Parses DI-tector output
	virDict = {}
	outVir = open(outVir, "r")
	
	line = outVir.readline() # Reads line
	line = outVir.readline() # Reads line
	lineList = line.split("\t") # Splits line into list

	for i in range(len(lineList[:-1])): # Final item in list is a "/n"
		item = lineList[i].split("_")
		bp = int(item[0])
		ri = int(item[2])

		virDict[(bp, ri)] = int(item[4])

	dataDict.update({"virDict": virDict})


def compile(outputPath, dataDict, sd): # puts into CSV alongside simulated reads

	output = open(outputPath, "w")
	
	simDict = dataDict["simDict"]
	ditDict = dataDict["ditDict"]
	virDict = dataDict["virDict"]

	totalDitCount = 0
	totalVirCount = 0

	unidList = []

	for key in simDict.keys():
		ditCount = 0
		virCount = 0

		filtDict= {k: v for k, v in ditDict.items() if (k[0]+sd >= key[0] and k[0]-sd <= key[0]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1])}
		for x in filtDict:
			ditCount += filtDict[x]
		
		filtDict= {k: v for k, v in virDict.items() if (k[0]+sd >= key[0] and k[0]-sd <= key[0]) and (k[1]+sd >= key[1] and k[1]-sd <= key[1])}
		for x in filtDict:
			virCount += filtDict[x]
		
		totalDitCount += ditCount
		totalVirCount += virCount

		output.write(str(key) + "\t" + str(simDict[key]) + "\t" + str(ditCount) + "\t" + str(virCount) + "\n")

		if ditCount == 0 and virCount == 0:
			unidList.append([key])

	output.close()


			

"""
Things to do:

> Summary file

> Be able to extract fasta sequence of unfound sequences

> Update readme

"""
			







if __name__ == "__main__":
	main(sys.argv[1:]) 

