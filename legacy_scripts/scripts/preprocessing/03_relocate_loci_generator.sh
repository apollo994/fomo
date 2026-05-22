for s in $(ls ../../input)
do
	gff=$(find ../../input/$s -name "${s}.lnc_RNA.longest.gff3")
	bed=$(find ../../input/$s -name "${s}.intergenic.bed")
	out_gff="${gff%.lnc_RNA.longest.gff3}.decoy.longest.gff3"
	echo python3 ./03_relocate_loci.py $(realpath $bed) $(realpath $gff) $(realpath -m $out_gff)
done
