for s in $(ls ../../input)
do
	fna=$(find ../../input/$s -name "*fna")
	gff=$(find ../../input/$s -name "*decoy.longest.gff3")
	out_fa="${gff%.gff3}.fa"
	echo "source ~/miniconda3/etc/profile.d/conda.sh && conda activate agat && agat_sp_extract_sequences.pl -g $(realpath $gff) -f $(realpath $fna) -t exon --merge -o $(realpath -m $out_fa)"
done
