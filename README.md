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

**Installation**:
<br>Ensure command is excecutable using chmod

<br>

---

**Usage**:
```
DIG.py <SimulationMethod> [options]
```
<br>

_Simulation Methods:_

Method | Description
:----- | :-----------
`MBP` | See below
`INDEL` | See below
`Copyback` | See below
`MultiSeg` | See below
`MultiSeg2` | See below
`NoDIP` | See below


<br>

_Compulsory Flags:_

Flag | Description
:--- | :-----------
`-f` / `--file` | Input fasta file (string)
`-o` / `--outdir` | Output directory (string)
`-m` / `--max` | Maximum read length (int)
`-t` / `--total` | Total number of reads (int)

<br>

_Method Specific Flags:_

Flag | Description
:--- | :-----------
`-c` / `--copybackratio` | Copyback only - sets ratio for 5' copyback, 5' snapback, 3' copyback, and 3' snapback DIPs (Default= t=0.45,0.05,0.45,0.05) 
`-n` / `--min` | MultiSeg only - sets minimum read length nucleotides (Defaul 300)
`--fragment` | MultiSeg only - Fragment reads, default = False
`-x` / `--num` | MultiSeg only - Number of fragments, default = 100000
`-l` / `--len` | MultiSeg only - Average read length, default = 300
`-s` / `--std` | MultiSeg only - Length standard deviation, default = 50)


<br>


---


**Methods**:

_MBP_
<br>Lengths with break point and reintination point in the middle of the resulting read with all reads at a given whole read length.
<br>Example: `DIG.py MBP -f reference.fasta -o outputDirectory -m 180 -t 100000`

_INDEL_
<br>Lengths with break point and reintination point in a random point in the resulting read with all reads at a given whole read length.
<br>Example: `DIG.py INDEL -f reference.fasta -o outputDirectory -m 180 -t 100000`

_Copyback_
<br>Generates copy back reads from 3' and 5', with a random segment in middle before copying back, also include snapback (no segment between reverse copied read) with all reads at a given whole read length.
<br>Example: `DIG.py Copyback -f reference.fasta -o outputDirectory -m 180 -t 100000 -c 0.45,0.05,0.45,0.05`

_MultiSeg_
<br>Reads created from first and last 600nt at random lengths within a given minimum and maximum, from random segments.
<br>Example with fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 1200 -t 1000 -n 300 -s 50 --fragment -x 100000 -l 300 -s 50`
<br>Example without fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 180 -t 100000 -n 300 -s 50`

_MultiSeg2_
<br>Reads created from first and last 600nt at random lengths within a given minimum and maximum, from random segments.
<br>Example: `DIG.py MultiSeg2 -f reference.fasta -o outputDirectory -m 150 -t 100000`

_NoDIP_
<br>Lengths created as just fragments of an existing genomes, i.e. reads with no dips.
<br>Example: `DIG.py NoDIP -f reference.fasta -o outputDirectory -m 150 -t 100000`


<br>

![DIG.py Diagrams](https://github.com/janewalls/DIG/blob/main/Images/DIG_diagrams.png?raw=true)

---

**Output**:

_Fasta files;_ <br> Saved in output directory e.g. SimMBP.fasta. Note: files will override if written in the same

_Summary files;_ <br>Output in csv saved in output directory. Items in csv saved as: read #, start, bp, ri, end

_Fragment Option;_ <br>(MultiSeg method only)Output in csv saved in output directory. Items in csv saved as: read #, fragment #, fragment length, bp boolean

<br>

---

**Output Parser Program:**

<br>

Program can be used to compare reads identified.

<br>

```
OutParse.py [options]
```

<br>

_Description:_

Used to compare outputs, has 2 usages, to compare outputs between simulations as created from DIG.py

<br>

_Flag Options:_

Flag | Description
:--- | :----------
`-s` / `--sim` | Simulated data output csv file (string)
`-d` / `--ditector` | Maximum read length (string)
`-v` / `--virema` | Total number of reads (string)
`-o` / `--outfile` | Output directory (string)
`--std` | Standard deviation for matches (int)
`-f` / `--fasta` | Simulated reads fasta file (string)

<br>

_Example:_

`OutParse.py -s simData -o outputDirectory -m 180 -t 100000 -n 300 -s 50`

<br>

---

**Example Bash Scripts:**



<br>

---

**Source**:

_MBP method;_ <br>Routh, A. and Johnson, J.E., 2014. Discovery of functional genomic motifs in viruses with ViReMa–a Virus Recombination Mapper–for analysis of next-generation sequencing data. Nucleic acids research, 42(2), pp.e11-e11.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902915/

_INDEL method_ & _Copyback method;_  <br>Beauclair, G., Mura, M., Combredet, C., Tangy, F., Jouvenet, N. and Komarova, A.V., 2018. DI-tector: defective interfering viral genomes’ detector for next-generation sequencing data. RNA, 24(10), pp.1285-1296.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140465/

_MultiSeg method;_ <br>Alnaji, F.G., Holmes, J.R., Rendon, G., Vera, J.C., Fields, C.J., Martin, B.E. and Brooke, C.B., 2019. Sequencing framework for the sensitive detection and precise mapping of defective interfering particle-associated deletions across influenza A and B viruses. Journal of virology, 93(11), pp.e00354-19.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6532088/

