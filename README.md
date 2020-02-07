# Still in development
# CNV_FunctionalAnnotation (All currently done in hg19)
Details to annotate CNV and apply prediction model

Required python3 dependencies: defaultdict, pybedtools, numpy, pandas, scipy

## Citation
Han L, Zhao X, Benton ML, Perumal T, Collins RL, Hoffman, GE, Johnson JS, Sloofman L, Stone MJ, CommonMind Consortium, Brennand KJ, Brand H, Sieberts SK, Marenco S, Peters MA, Lipska BK, Roussos P, Capra JA, Talkowski M, Ruderfer DM. Functional annotation of rare structural variation in the human brain. BioRxiv doi: https://doi.org/10.1101/711754

## Data Generation

To download/generate annotation files run `bin/make_annotation_files`.

## Annotating CNV

perl bin/annotate-cnv.pl [input_file] > [output_file]

There is a test CNV file in test/test.cnv (see below). Only DEL or DUP are accepted types
id chr start end type  
CNV1 2 133450000 133550000 DEL  
CNV2 12 12150000 12880000 DEL  
CNV3 7 6450000 6550000 DUP  
CNV4 4 70450000 71260000 DUP  
CNV5 18 5450000 7550000 DEL  

Output will include all genes affected by any CNV as separate rows. CNVs that do not affect genes will not be output currently.   
Output file has same columns as input file plus below:  
len: CNV length  
gene: Ensembl gene name  
exon: proportion of exonic sequence affected   
pro: proportion of promoter sequence affected  
enh: sum proportion of all enhancers sequence affected   
intad: within TAD (1) or not (0)  

## Prediction models

The deletion and duplication models are saved in CNV-models.Rdata


## Calculating regulatory disruption score
