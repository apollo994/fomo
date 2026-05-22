INPUT_DIR="../../benchmark/input"
GFFCOMPARE_DIR="../../benchmark/gffcompare"
OUT_DIR="../../benchmark/TPFP_analysis/input"

for gff in $(find "$INPUT_DIR" -name "*.gff3" -not -name "*_decoy.gff3")
do
	base=$(basename "$gff" .gff3)
	annot="$GFFCOMPARE_DIR/$base/$base.annotated.gtf"
	tp="$OUT_DIR/${base}.TP.gff3"
	fp="$OUT_DIR/${base}.FP.gff3"
	echo "python3 ./05_split_TPFP.py $(realpath $gff) $(realpath $annot) $(realpath -m $tp) $(realpath -m $fp)"
done
