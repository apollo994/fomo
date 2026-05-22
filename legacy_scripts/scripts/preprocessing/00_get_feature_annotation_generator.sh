for s in $(ls ../../input)
do
	fna=$(find ../../input/$s -name "*fna")
	gff=$(find ../../input/$s -name "*gff3")
	echo bash ./00_get_feature_annotation.sh $(realpath $gff) $(realpath $fna) mRNA
	echo bash ./00_get_feature_annotation.sh $(realpath $gff) $(realpath $fna) lnc_RNA
done
