ref_dir="../../input/Polyommatus_icarus_265386/GCA_937595015.1"
out_root="../../benchmark/gffcompare"

for gff in $(find ../../benchmark/input -name "*.gff3")
do
	base=$(basename "$gff" .gff3)
	genetype=${base##*_}
	case $genetype in
		pcg)   ref="$ref_dir/Polyommatus_icarus_265386.mRNA.longest.gff3" ;;
		lnc)   ref="$ref_dir/Polyommatus_icarus_265386.lnc_RNA.longest.nonoverlapping.gff3" ;;
		decoy) ref="$ref_dir/Polyommatus_icarus_265386.decoy.longest.gff3" ;;
		*)     echo "skip: unknown genetype in $base" >&2; continue ;;
	esac
	outdir="$out_root/$base"
	echo "mkdir -p $(realpath -m $outdir) && gffcompare -M --no-exon-merge -r $(realpath $ref) -o $(realpath -m $outdir)/$base $(realpath $gff)"
done
