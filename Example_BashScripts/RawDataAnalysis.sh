# Raw Data Analysis & Quality Control


# Put all raw read files in a list 
sampleList=""
for sample in Bld_33_S1_R1_001 Bld_33_S1_R2_001  Bld_36_S2_R1_001 Bld_36_S2_R2_001 Bld_38_S4_R1_001 Bld_38_S4_R2_001 Bld_40_S6_R1_001 Bld_40_S6_R2_001 Isol_19_S17_R1_001 Isol_19_S17_R2_001
do
	sampleList="${sampleList} ${data}/${sample}.fastq.bz2"
done

fastqc -o $fastqc_dir ${sampleList}  # Quality control all reads in list


for sample in Bld_33_S1 Bld_36_S2 Bld_38_S4 Bld_40_S6 Isol_19_S17 
do
	r1="${data}/${sample}_1.fastq.bz2"
	r2="${data}/${sample}_2.fastq.bz2"

	trim_galore --paired ${r1} ${r2} -o ${trimdata} # Quality trim data

	trim1="${trimdata}/${sample}_1.fastq.bz2_val_1.fq"
	trim2="${trimdata}/${sample}_2.fastq.bz2_val_2.fq"
	samout="${outdir}/${sample}.sam"
	bamsorted="${outdir}/${sample}.sorted.bam"	

	bowtie2 --phred33 -x $bowtieIndex -1 ${trim1} -2 ${trim2} -S ${samout} 2> ${sample}.log  # Align to BTV

	samtools view -buS ${samout} | samtools sort -o ${bamsorted} 

	samtools index ${bamsorted}

	getinsertsize.py ${samout}
done


# Align to cow data

sampleList=""

bowtie2-build $reference $bowtieIndex # Create cow bowtie2 index

for sample in Bld_33_S1 Bld_36_S2 Bld_38_S4 Bld_40_S6 Isol_19_S17
do
	trim1="${trimdata}/${sample}_1.fastq.bz2_val_1.fq"
	trim2="${trimdata}/${sample}_2.fastq.bz2_val_2.fq"
	samout="${outdir}/${sample}.sam"
	bamsorted="${outdir}/${sample}.sorted.bam"	

	bowtie2 --phred33 -x $bowtieIndex -1 ${trim1} -2 ${trim2} -S ${samout} 2> ${sample}.log

	samtools view -buS ${samout} | samtools sort -o ${bamsorted} 

	samtools index ${bamsorted}

	getinsertsize.py ${samout} > ${sample}_insert.log

	bamCoverage -b $bamsorted -of bedgraph -o ${outdir}/stats/${sample}.bed
	bamCoverage -b $bamsorted -of bigwig -o ${outdir}/stats/${sample}.bw


	input="${data}/${sample}.sam"
	output="${outdir}/${sample}.unmapped.fastq"

	# Filter out cow data
	samtools fastq -f 4 $input > $output

	sampleList="${sampleList} ${data}/${sample}.unmapped.fastq" # put filtered sample in list
done

# QC on filtered reads

fastqc -o $fastqc_dir ${sampleList}










