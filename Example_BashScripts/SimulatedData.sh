# Simulate and Analysis Data

# Simulate data
DIG.py MBP -f $reference -o $data -m 150 -t 100000 --seg 1

DIG.py INDEL -f $reference -o $data -m 150 -t 100000 --seg 10

DIG.py CopyBack -f $reference -o $data -m 150 -t 100000 -c 0.45,0.05,0.45,0.05 --seg 8

DIG.py MultiSeg2 -f $reference -o $data -m 150 -t 100000 -w 600

DIG.py NoDIP -f $reference -o $data -m 150 -t 100000 --seg 5


# Create /usr/bin/time output files
to_dit="${output}/timeOutput_DIT.txt"
to_vir="${output}/timeOutput_VRM.txt"
> $to_dit
> $to_vir


# Analyse data
for samplef in SimCopyBack SimIndel SimMBP SimMultiSeg2 SimNoDIP
do
	path1="${data}/${samplef}.fasta"
	sim="${data}/${samplef}.csv"
	ditecor="${output}/${samplef}/${samplef}_DIT/DI-tector_output_sorted.txt"
	virema="${output}/${samplef}/${samplef}_VRM/Virus_Recombination_Results.txt"
	outdir="${output}/${samplef}"

	# Create directories
	mkdir -p $outdir
	mkdir -p ${outdir}/res_std

	echo "${samplef}\n" >> $to_dit
	echo "${samplef}\n" >> $to_vir

	/usr/bin/time -f "time\t%e\navgMem\t%K\ncpu\t%P\nmaxRes\t%M\n" -ao $to_vir ViReMa.py $reference $path1 ViReMa_sim4_${samplef}.txt --p 4 --Aligner bwa -F --Output_Dir ${outdir}/${samplef}_VRM 
	
	/usr/bin/time -f "time\t%e\navgMem\t%K\ncpu\t%P\nmaxRes\t%M\n" -ao $to_dit python3 DI-tector_06.py $reference $path1 -fk -n 1 -x 4 -o ${outdir}/${samplef}_DIT

	OutParse.py -s $sim -d $ditecor -v $virema -o $outdir 
	OutParse.py -s $sim -d $ditecor -v $virema -o ${outdir}/res_std --std 2
done







