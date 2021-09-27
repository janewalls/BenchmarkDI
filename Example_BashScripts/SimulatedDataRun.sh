
DIG.py MBP -f $reference -o $data -m 150 -t 100000 --seg 1

DIG.py INDEL -f $reference -o $data -m 150 -t 100000 --seg 10

DIG.py CopyBack -f $reference -o $data -m 150 -t 100000 -c 0.45,0.05,0.45,0.05 --seg 8

DIG.py MultiSeg2 -f $reference -o $data -m 150 -t 100000 -w 600

DIG.py NoDIP -f $reference -o $data -m 150 -t 100000 --seg 5s


for sample in SimCopyBack SimIndel SimMBP SimMultiSeg2 SimNoDIP
do
	mkdir -p $outdir
	mkdir -p ${outdir}/res_std

	ViReMa.py $reference $path ViReMa_sim3_${sample}.txt --Aligner bwa -F --Output_Dir ${outdir}/${sample}_VRM 
	
	python3 $DITector $reference $path -fk -n 1 -x 3 -o ${outdir}/${sample}_DIT

	OutParse.py -s $sim -d $ditecor -v $virema -o $outdir
	OutParse.py -s $sim -d $ditecor -v $virema -o ${outdir}/res_std --std 2
done
