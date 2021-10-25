# BTV DIP Detection Analysis


# Analyse filtered 
for sample in Bld_33_S1 Bld_36_S2 Bld_38_S4 Bld_40_S6 Isol_19_S17
do
	path1="${data}/${sample}.unmapped.fastq"
	outdir="${output}/${sample}"
	ditecor="${outdir}/${sample}_DIT/DI-tector_output_sorted.txt"
	virema="${outdir}/${sample}_VRM/Virus_Recombination_Results.txt"

	mkdir -p $outdir
	mkdir -p ${outdir}/res_std

	ViReMa.py $reference $path1 ViReMa_${sample}.txt --Aligner bwa --bed --Output_Dir ${outdir}/${sample}_VRM 
	
	python3 DI-tector_06.py $reference $path1 -k -n 1 -x 3 -o ${outdir}/${sample}_DIT

	OutParse.py -d $ditecor -v $virema -o ${outdir}/res_std --std 1 # Compare tool outputs
	
	ditToBed.py -d ${outdir}/${sample}_DIT # Create bed file for DI-Tecor
done

# Compile sample outputs
CompDIP.py --s1 ${output}/Bld_33_S1/res_std/parser_output.csv --s2 ${output}/Bld_36_S2/res_std/parser_output.csv --s3 ${output}/Bld_38_S4/res_std/parser_output.csv --s4 ${output}/Bld_40_S6/res_std/parser_output.csv -o ${output}/compiledOutput.txt
CompDIP.py --s1 ${output}/Bld_33_S1/res_std/parser_output.csv --s2 ${output}/Bld_36_S2/res_std/parser_output.csv --s3 ${output}/Bld_38_S4/res_std/parser_output.csv --s4 ${output}/Bld_40_S6/res_std/parser_output.csv -o ${output}/compiledOutput_cut.txt --cut 1








