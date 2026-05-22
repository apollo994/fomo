outdir="../../benchmark/summary_stats"
for gff in $(find ../../benchmark/input -name "*.gff3")
do
	base=$(basename "$gff" .gff3)
	without_tool=${base#*_}
	target=${without_tool%%_from_*}
	fna=$(find ../../input/$target -name "*fna")
	out="$outdir/${base}.stats.txt"
	yaml="$outdir/${base}.stats.yaml"
	echo "source ~/miniconda3/etc/profile.d/conda.sh && conda activate agat && agat_sp_statistics.pl -i $(realpath $gff) -g $(realpath $fna) -o $(realpath -m $out) --yaml $(realpath -m $yaml)"
done
