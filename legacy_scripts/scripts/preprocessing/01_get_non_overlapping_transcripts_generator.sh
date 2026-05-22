for s in $(ls ../../input)
do
	fna=$(find ../../input/$s -name "*fna")
	gff1=$(find ../../input/$s -name "*lnc_RNA.longest.gff3")
	gff2=$(find ../../input/$s -name "*mRNA.longest.gff3")
	echo bash ./01_get_non_overlapping_transcripts.sh $(realpath $gff1) $(realpath $gff2) $(realpath $fna)
done
