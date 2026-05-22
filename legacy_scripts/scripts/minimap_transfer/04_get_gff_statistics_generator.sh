for d in $(ls ../../minimap_models)
do
	target=${d%_from_*}
	fna=$(find ../../input/$target -name "*fna")
	for gff in $(find ../../minimap_models/$d -name "*.gff3")
	do
		out="${gff%.gff3}.stats.txt"
		yaml="${gff%.gff3}.stats.yaml"
		echo "source ~/miniconda3/etc/profile.d/conda.sh && conda activate agat && agat_sp_statistics.pl -i $(realpath $gff) -g $(realpath $fna) -o $(realpath -m $out) --yaml $(realpath -m $yaml)"
	done
done
