# DIG
### -> **DI G**enerator

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
python3 DIG.py SimulationMethod [options]
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
`-f` / `--file` | Input fasta file
`-o` / `--outdir` | Output directory
`-m` / `--max` | Maximum read length
`-t` / `--total` | Total number of reads

<br>

_Method Specific Flags:_
Flag | Description
-------|------------
`-cbr` / `--copybackratio` | Copyback only - sets ratio for 5' copyback, 5' snapback, 3' copyback, and 3' snapback DIPs (Default= 0.45,0.05,0.45,0.05)
`-n` / `--min` | MultiSeg only - sets minimum read length (Default= 300 nt)


<br>


---


**Methods**:

_ViReMa_
<br>Lengths with break point and reintination point 

_INDEL_
<br>Lengths with 

_Copyback_
<br>Generates copy back reads from 3' and 5', can also include snapback

_MultiSeg_
<br>Fragments occur

<br>

---

**Source**:

_ViReMa method;_ <br>Routh, A. and Johnson, J.E., 2014. Discovery of functional genomic motifs in viruses with ViReMa–a Virus Recombination Mapper–for analysis of next-generation sequencing data. Nucleic acids research, 42(2), pp.e11-e11.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3902915/

_INDEL method_ & _Copyback method;_  <br>Beauclair, G., Mura, M., Combredet, C., Tangy, F., Jouvenet, N. and Komarova, A.V., 2018. DI-tector: defective interfering viral genomes’ detector for next-generation sequencing data. RNA, 24(10), pp.1285-1296.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140465/

_MultiSeg method_ - <br>Alnaji, F.G., Holmes, J.R., Rendon, G., Vera, J.C., Fields, C.J., Martin, B.E. and Brooke, C.B., 2019. Sequencing framework for the sensitive detection and precise mapping of defective interfering particle-associated deletions across influenza A and B viruses. Journal of virology, 93(11), pp.e00354-19.<br>https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6532088/

