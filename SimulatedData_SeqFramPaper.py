import random
import Bio

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Adjustable parameters 
fastaFileS1 = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
fastaFileS2 = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/example.fasta"
maxLength = 188
totalReads = 50
count = 0

readOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.fasta"
summaryOutput = "/Users/janewalls/Documents/VS_CODE/MastersProjectDI/SimDITecINDEL.csv"
