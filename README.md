1;95;0c# Still in development
# CNV_FunctionalAnnotation (All currently done in hg19)
Details to annotate CNV and apply prediction model

Required python3 dependencies: defaultdict, pybedtools, numpy, pandas, scipy

## Citation
Han L, Zhao X, Benton ML, Perumal T, Collins RL, Hoffman, GE, Johnson JS, Sloofman L, Stone MJ, CommonMind Consortium, Brennand KJ, Brand H, Sieberts SK, Marenco S, Peters MA, Lipska BK, Roussos P, Capra JA, Talkowski M, Ruderfer DM. Functional annotation of rare structural variation in the human brain. BioRxiv doi: https://doi.org/10.1101/711754

## Data Generation

To download/generate annotation files run `bin/make_annotation_files`.

## Prediction models

The deletion and duplication models are saved in CNV-models.Rdata


## Annotating CNV and applying models to generate regulatory disruption scores (will generate intermediate and results files in same directory as CNV file)

sh bin/annotate-cnv.sh [path where git repository was cloned] [path_to_CNV_file]

There is a test CNV file in test/test.cnv (see below). Only DEL or DUP are accepted types

id  |chr| start| end|SVType  
---|---|---|---|---
CNV1|2|133450000|133550000| DEL  
CNV2|12|12150000|12880000|DEL   
CNV3|7|6450000|6550000|DUP  
CNV4|4|70450000|71260000|DUP   
CNV5|18|5450000|7550000|DEL   


## Intermediate CNV by gene annotation file

Output file has same columns as input file plus below:  

*SVLen*: CNV length  
*gene*: Ensembl gene name  
*trans*: proportion of entire gene affected
*ExonProp*: proportion of exonic sequence affected   
*pro.sum*: proportion of promoter sequence affected  
*enh.sum*: sum proportion of all enhancers sequence affected   
*inTad*: within TAD (1) or not (0)  

id |chr |start |end |SVType |SVLen |gene |trans |ExonProp |pro.sum |enh.sum |inTad
---|---|---|---|---|---|---|---|---|---|---|---
CNV1 |2 |133450000 |133550000 |DEL |100000 |ENSG00000176771 |0.111525 |0.50330015715034 |0 |0 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000070018 |1 |1 |1 |0 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000111261 |1 |1 |1 |3 |0
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000111266 |1 |1 |1 |0 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000111269 |1 |1 |1 |0 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000111276 |1 |1 |1 |6 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000121380 |1 |1 |1 |5 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000165714 |1 |1 |1 |0 |1
CNV2 |12 |12150000 |12880000 |DEL |730000 |ENSG00000178878 |0.0110419 |0.0354825348612453 |1 |0 |1


## Final output file

id |chr |start |end |SVType |SVLen |ngenes |exon |enh |pro |pred_exp |reg_dist
---|---|---|---|---|---|---|---|---|---|---|---
CNV1 |2 |133450000 |133550000 |DEL |100000 |1 |0.50330015715034 |0 |0 |-0.923467521409079 |-0.683216839007316
CNV2 |12 |12150000 |12880000 |DEL |730000 |22 |21.0354825348612 |36 |22 |-45.6165484017837 |-14.6298854256054
CNV3 |7 |6450000 |6550000 |DUP |100000 |5 |4.34803581412534 |21 |4 |4.90544957250558 |1.67441995009775
CNV4 |4 |70450000 |71260000 |DUP |810000 |21 |20.4703337453646 |0 |21 |30.8031235891722 |9.62781209180071
CNV5 |18 |5450000 |7550000 |DEL |2100000 |37 |36.2654214123007 |2 |37 |-88.1485171766695 |-8.63831141801943

