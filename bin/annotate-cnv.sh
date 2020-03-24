
input=$1

echo $1


echo "creating CNV bed file (assuming first line is header)..."
awk 'NR>1{OFS=""; print "chr",$2,"\t",$3,"\t",$4,"\t",$1}' $input > ${input}.bed

echo "Characterizing genic and exonic effects..."
bedtools intersect -a ${input}.bed -b Ensemblv75/CMC-gene.bed -wo | awk -F '\t' '$7-$6 > 0{OFS="\t"; print $4,$8,$10/($7-$6)}' > ${input}-gene
bedtools intersect -a ${input}.bed -b Ensemblv75/CMC-gene-exon.bed -wo | awk -F '\t' '{OFS="\t"; print $4,$8,$9+1}' > ${input}-exon

echo "Characterizing TAD effects..."
bedtools intersect -a ${input}.bed -b chromatin_structure/DER-18_TAD_adultbrain.bed -wo | awk -F '\t' '{OFS=""; print $4,"\t",$5,":",$6,"-",$7}' > ${input}-tad
bedtools intersect -a Ensemblv75/CMC-gene.bed -b chromatin_structure/DER-18_TAD_adultbrain.bed -wo | awk -F '\t' '{OFS=""; print $4,"\t",$6,":",$7,"-",$8}' > ${input}-gene-tad

echo "Characterizing enhancer effects..."
bedtools intersect -a ${input}.bed -b enhancers/DER-03a_hg19_PEC_enhancers.bed -wo | awk '{OFS="\t"; print $0,$9/($7-$6)}' > ${input}-PEC-enhancers.all
awk '$10 == 1 {print $4}' ${input}-PEC-enhancers.all | sort | uniq -c | awk '{print $2,$1}' > ${input}-PEC-enhancers.cnt
bedtools intersect -a ${input}.bed -b enh-gene_links/PEC-enh-gene-links.bed -wo | awk -F '\t' '{OFS="\t"; print $4,$9,$10/($7-$6)}' > ${input}-PEC-enhancer-predictions

echo "Characterizing CTCF effects..."
bedtools intersect -a ${input}.bed -b ctcf/brain_tissues/ctcf_brain_merged.bed -wa -wb | awk '{print $4}' | sort | uniq -c | awk '{print $2,$1}' > ${input}-ctcf.cnt

echo "Characterizing promoter effects..."
awk -F '\t' '{OFS="\t";if($5 == "+" && $2 > 2000){print $1,$2-2000,$2,$4};if($5 == "+" && $2 < 2000){print $1,0,$2,$4};if($5 == "+"){print $1,$3,$3+2000,$4}}' Ensemblv75/CMC-gene.bed > ${input}-gene-promoters
bedtools intersect -a ${input}.bed -b promoters/CMC-gene-promoters -wo | awk -F '\t' '{OFS="\t"; print $4,$8,$9/($7-$6)}' > ${input}-promoters


