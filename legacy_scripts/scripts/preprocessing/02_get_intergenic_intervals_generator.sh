for s in $(ls ../../input)
do
	gff=$(find ../../input/$s -name "${s}.gff3")
	out_bed="${gff%.gff3}.intergenic.bed"
	echo bash ./02_get_intergenic_intervals.sh $(realpath $gff) $(realpath -m $out_bed)
done
