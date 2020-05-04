
dir=$1
input=$2

echo "CNV_Functional_Annotation Directory: " $1
echo "CNV file with path: " $2 


echo "creating CNV bed file (assuming first line is header)..."
awk 'NR>1{OFS=""; print "chr",$2,"\t",$3,"\t",$4,"\t",$1}' $input > ${input}.bed

echo "Characterizing genic and exonic effects..."
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/Ensemblv75/CMC-gene.bed -wo | awk -F '\t' '$7-$6 > 0{OFS="\t"; print $4,$8,$10/($7-$6)}' > ${input}-gene
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/Ensemblv75/CMC-gene-exon.bed -wo | awk -F '\t' '{OFS="\t"; print $4,$8,$9+1}' > ${input}-exon

echo "Characterizing TAD effects..."
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/chromatin_structure/DER-18_TAD_adultbrain.bed -wo | awk -F '\t' '{OFS=""; print $4,"\t",$5,":",$6,"-",$7}' > ${input}-tad
bedtools intersect -a ${dir}/CNV_FunctionalAnnotation/Ensemblv75/CMC-gene.bed -b ${dir}/CNV_FunctionalAnnotation/chromatin_structure/DER-18_TAD_adultbrain.bed -wo | awk -F '\t' '{OFS=""; print $4,"\t",$6,":",$7,"-",$8}' > ${input}-gene-tad

echo "Characterizing enhancer effects..."
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/enhancers/DER-03a_hg19_PEC_enhancers.bed -wo | awk '{OFS="\t"; print $0,$9/($7-$6)}' > ${input}-PEC-enhancers.all
awk '$10 == 1 {print $4}' ${input}-PEC-enhancers.all | sort | uniq -c | awk '{print $2,$1}' > ${input}-PEC-enhancers.cnt
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/enh-gene_links/PEC-enh-gene-links.bed -wo | awk -F '\t' '{OFS="\t"; print $4,$9,$10/($7-$6)}' > ${input}-PEC-enhancer-predictions

echo "Characterizing CTCF effects..."
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/ctcf/brain_tissues/ctcf_brain_merged.bed -wa -wb | awk '{print $4}' | sort | uniq -c | awk '{print $2,$1}' > ${input}-ctcf.cnt

echo "Characterizing promoter effects..."
awk -F '\t' '{OFS="\t";if($5 == "+" && $2 > 2000){print $1,$2-2000,$2,$4};if($5 == "+" && $2 < 2000){print $1,0,$2,$4};if($5 == "+"){print $1,$3,$3+2000,$4}}' ${dir}/CNV_FunctionalAnnotation/Ensemblv75/CMC-gene.bed > ${input}-gene-promoters
bedtools intersect -a ${input}.bed -b ${dir}/CNV_FunctionalAnnotation/promoters/Ensemblv75-promoters -wo | awk -F '\t' '{OFS="\t"; print $4,$8,$9/($7-$6)}' > ${input}-promoters

echo "Running R script to apply model..."
Rscript ${dir}/CNV_FunctionalAnnotation/bin/apply-model.R $1 $2

echo "Deleting intermediate files..."
rm ${input}-gene ${input}-exon ${input}-tad ${input}-gene-tad ${input}-PEC-enhancers.all ${input}-PEC-enhancers.cnt ${input}-PEC-enhancer-predictions ${input}-ctcf.cnt ${input}-gene-promoters ${input}-promoters 


