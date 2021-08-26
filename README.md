# DIG
### :dna: *DI G*enerator :petri_dish:

A program to generate Defective Interfering particles with 4 different methods


<br>

---

**Requirements**:
<br>Python 3.0+ 
<br>BioPython

<br>

---

**Usage**:
```
DIG.py <SimulationMethod> [options]
```
<br>

_Simulation Methods:_
Method | Description
-------|------------
`ViReMa` | See below
`INDEL` | See below
`Copyback` | See below
`MultiSeg` | See below

<br>

_Compulsory Flags:_
Flag | Description
-------|------------
`-f` / `--file` | Input fasta file (string)
`-o` / `--outdir` | Output directory (string)
`-m` / `--max` | Maximum read length (int)
`-t` / `--total` | Total number of reads (int)

<br>

_Method Specific Flags:_
Flag | Description
-------|------------
`-c` / `--copybackratio` | Copyback only - sets ratio for 5' copyback, 5' snapback, 3' copyback, and 3' snapback DIPs (Default= 0.45,0.05,0.45,0.05) 
`-n` / `--min` | MultiSeg only - sets minimum read length nucleotides (Default= 300)
`--fragment` | MultiSeg only - Fragment reads, default = False
`-x` / `--num` | MultiSeg only - Number of fragments, default = 100000
`-l` / `--len` | MultiSeg only - Average read length, default = 300
`-s` / `--std` | MultiSeg only - Length standard deviation, default = 50)


<br>


---


**Methods**:

_ViReMa_
<br>Lengths with break point and reintination point in the middle of the resulting read with all reads at a given whole read length.
<br>Example: `DIG.py ViReMa -f reference.fasta -o outputDirectory -m 180 -t 50`

_INDEL_
<br>Lengths with break point and reintination point in a random point in the resulting read with all reads at a given whole read length.
<br>Example: `DIG.py INDEL -f reference.fasta -o outputDirectory -m 180 -t 50`

_Copyback_
<br>Generates copy back reads from 3' and 5', with a random segment in middle before copying back, also include snapback (no segment between reverse copied read) with all reads at a given whole read length.
<br>Example: `DIG.py Copyback -f reference.fasta -o outputDirectory -m 180 -t 50`

_MultiSeg_
<br>Reads created from first and last 600nt at random lengths within a given minimum and maximum, from random segments.
<br>Example with fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 180 -t 50`
<br>Example without fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 180 -t 50`

<br>

---

**Output**:

_Fasta files;_ <br>

_Summary files;_ <br>Output in csv saved in output directory. Items saved as 

<br>

---

**Source**:

_ViReMa method;_ <br>Routh, A. and Johnson, J.E., 2014. Discovery of functional genomic motifs in viruses with ViReMa–a Virus Recombination Mapper–for analysis of next-generation sequencing data. Nucleic acids research, 42(2), pp.e11-e11.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902915/

_INDEL method_ & _Copyback method;_  <br>Beauclair, G., Mura, M., Combredet, C., Tangy, F., Jouvenet, N. and Komarova, A.V., 2018. DI-tector: defective interfering viral genomes’ detector for next-generation sequencing data. RNA, 24(10), pp.1285-1296.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140465/

_MultiSeg method;_ <br>Alnaji, F.G., Holmes, J.R., Rendon, G., Vera, J.C., Fields, C.J., Martin, B.E. and Brooke, C.B., 2019. Sequencing framework for the sensitive detection and precise mapping of defective interfering particle-associated deletions across influenza A and B viruses. Journal of virology, 93(11), pp.e00354-19.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6532088/

