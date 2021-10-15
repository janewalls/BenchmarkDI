# DIG
### :dna: *DI G*enerator :petri_dish:

A program to generate Defective Interfering particles with 4 different methods, and other packages to help with analysis the output.
<br>

There were 4 programs created:
- DIG = Simulate Defective Particles
- OutParse = Compare DIPs found between tools, and if present simulated DIPs
- CompDIP = compile and compare DIPs in samples, and tools
- ditToBed = Create bed file from DI-Tector output


<br>

## Table of Contents:

- [Requirements](#Requirements)
- [Installation](#Installation)
- [DIG.py](#DIG.py)
    - Usage
    - Methods
    - Output
    - Example
- [Outparse.py](#Outparse.py)
    - Usage
    - 
- [CompDIP.py](#CompDIP.py)
    - Usage
- [DitToBed.py](#DitToBed.py)
    - Usage

<br>

---

### Requirements:
<br>Python 3.0+ 
<br>BioPython

<br>

---

### Installation:
<br>
Using the command line:

```
git clone https://github.com/janewalls/DIG.git
```

Then ensure all commands are executable: 

```
chmod 555 DIG.py OutParse.py CompDIP.py ditToBed.py
```
<br>
<br>
<br>

## DIG.py 

<br>

Simulate defective particles 

<br>

### Usage:
```
DIG.py <SimulationMethod> [options]
```
<br>

*_Simulation Methods:_*

Method | Description
:----- | :-----------
`MBP` | See below
`INDEL` | See below
`Copyback` | See below
`MultiSeg` | See below
`MultiSeg2` | See below
`NoDIP` | See below


<br>

*_Compulsory Flags:_*

Flag | Description | Type | Required
:--- | :---------- | :--- | --------:
`-f` / `--file` | Input fasta file | `String` | __Yes__
`-o` / `--outdir` | Output directory | `String` | __Yes__
`-m` / `--max` | Maximum read length | `int` | __Yes__
`-t` / `--total` | Total number of reads | `int` | __Yes__

<br>

*_Method Specific Flags:_*

Flag | Description | Type | Required
:--- | :---------- | :--- | --------:
`-c` / `--copybackratio` | Copyback only - sets ratio for 5' copyback, 5' snapback, 3' copyback, and 3' snapback DIPs (Default= 0.45,0.05,0.45,0.05)  | `String` | No
`-n` / `--min` | MultiSeg only - sets minimum read length nucleotides (Defaul 300) | `int` | No
`--fragment` | MultiSeg only - Fragment reads, default = False  | - | No
`-x` / `--num` | MultiSeg only - Number of fragments, default = 100000 | `int` | No
`-l` / `--len` | MultiSeg only - Average read length, default = 300 | `int` | No
`-s` / `--std` | MultiSeg only - Length standard deviation, default = 50) | `int` | No


<br>


---
<br>

### Methods:
<br>

_MBP_
<br>Lengths with break point and reintination point in the middle of the resulting read with all reads at a given whole read length.


_INDEL_
<br>Lengths with break point and reintination point in a random point in the resulting read with all reads at a given whole read length.


_Copyback_
<br>Generates copy back reads from 3' and 5', with a random segment in middle before copying back, also include snapback (no segment between reverse copied read) with all reads at a given whole read length.


_MultiSeg_
<br>Reads created from first and last 600nt at random lengths within a given minimum and maximum, from random segments.


_MultiSeg2_
<br>Reads created from first and last 600nt at random lengths within a given minimum and maximum, from random segments.


_NoDIP_
<br>Lengths created as just fragments of an existing genomes, i.e. reads with no dips.

<br>

![DIG.py Diagrams](https://github.com/janewalls/DIG/blob/main/Images/DIG_diagrams.png?raw=true)

---
<br>

### Examples: 

_MBP_
<br>`DIG.py MBP -f reference.fasta -o outputDirectory -m 180 -t 100000`

_INDEL_
<br>`DIG.py INDEL -f reference.fasta -o outputDirectory -m 180 -t 100000`

_Copyback_
<br>`DIG.py Copyback -f reference.fasta -o outputDirectory -m 180 -t 100000 -c 0.45,0.05,0.45,0.05`

_MultiSeg_
<br>With fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 1200 -t 1000 -n 300 -s 50 --fragment -x 100000 -l 300 -s 50`
<br>Without fragmentation: `DIG.py MultiSeg -f reference.fasta -o outputDirectory -m 180 -t 100000 -n 300 -s 50`

_MultiSeg2_
<br>`DIG.py MultiSeg2 -f reference.fasta -o outputDirectory -m 150 -t 100000`

_NoDIP_
<br>`DIG.py NoDIP -f reference.fasta -o outputDirectory -m 150 -t 100000`


---
<br>

### Output:
<br>

_Fasta files;_ <br> Saved in output directory e.g. SimMBP.fasta. Note: files will override if written in the same

_Summary files;_ <br>Output in csv saved in output directory. Items in csv saved as: read #, start, bp, ri, end

_Fragment Option;_ <br>(MultiSeg method only)Output in csv saved in output directory. Items in csv saved as: read #, fragment #, fragment length, bp boolean

<br>
<br>
<br>

## Outparse.py

<br>

Used to compare outputs, has 2 usages, to compare outputs between simulations as created from DIG.py

<br>

### Usage:

```
OutParse.py [options]
```

<br>



*_Flag Options:_*

Flag | Description | Type | Required
:--- | :---------- | :--- | --------:
`-s` / `--sim` | Simulated data output csv file | `String` | No
`-d` / `--ditector` | Maximum read length | `String` | __Yes__
`-v` / `--virema` | Total number of reads | `String` | __Yes__
`-o` / `--outfile` | Output directory | `String` | __Yes__
`--std` | Standard deviation for matches | `String` | No
`-f` / `--fasta` | Simulated reads fasta file | `String` | No

<br>

### Example:

`OutParse.py -s simData -o outputDirectory -m 180 -t 100000 -n 300 -s 50`

<br>
<br>
<br>

## CompDIP.py

<br>

Program to compile up to 5 sample results from OutParse.py, designed to only compare real samples, rather than simultation.
Minimum of 2 samples, maximum of 5

<br>

### Usage:

```
CompDIP.py [options]
```
<br>

_Flag Options:_

Flag | Description | Type | Required
:--- | :---------- | :--- | --------:
`--s1` | Sample 1 parser_out.csv file | `String` | __Yes__
`--s2` | Sample 2 parser_out.csv file | `String` | __Yes__
`--s3` | Sample 3 parser_out.csv file | `String` | No
`--s4` | Sample 4 parser_out.csv file | `String` | No
`--s5` | Sample 5 parser_out.csv file | `String` | No
`-o` | Output directory | `String` | __Yes__
`--cut` | Simulated reads fasta file | `int` | No

<br>
<br>
<br>

## ditToBed.py

<br>

Program to

<br>

### Usage:

```
ditToBed.py -d <Directory>
```
<br>

_Flag Options:_

Flag | Description | Type | Required
:--- | :---------- | :--- | --------:
`-d` | Directory of DI-Tector output | `String` | __Yes__


<br>

---

### Example Bash Scripts:
<br>

See [Example_BashScripts](/janewalls/DIG/Example_BashScripts) folder in GitHub

<br>

---

### References:

_(1) MBP method;_ <br>Routh, A. and Johnson, J.E., 2014. Discovery of functional genomic motifs in viruses with ViReMa–a Virus Recombination Mapper–for analysis of next-generation sequencing data. Nucleic acids research, 42(2), pp.e11-e11.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902915/

_(2) INDEL method_ & _Copyback method;_  <br>Beauclair, G., Mura, M., Combredet, C., Tangy, F., Jouvenet, N. and Komarova, A.V., 2018. DI-tector: defective interfering viral genomes’ detector for next-generation sequencing data. RNA, 24(10), pp.1285-1296.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140465/

_(3) MultiSeg method;_ <br>Alnaji, F.G., Holmes, J.R., Rendon, G., Vera, J.C., Fields, C.J., Martin, B.E. and Brooke, C.B., 2019. Sequencing framework for the sensitive detection and precise mapping of defective interfering particle-associated deletions across influenza A and B viruses. Journal of virology, 93(11), pp.e00354-19.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6532088/

