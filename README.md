# CNV_FunctionalAnnotation (All currently done in hg19)
Details to annotate CNV and apply prediction model

Required python3 dependencies: defaultdict, pybedtools, numpy, pandas, scipy

## Citation
Han L, Zhao X, Benton ML, Perumal T, Collins RL, Hoffman, GE, Johnson JS, Sloofman L, CommonMind Consortium, Brennand KJ, Brand H, Sieberts SK, Marenco S, Peters MA, Lipska BK, Roussos P, Capra JA, Talkowski M, Ruderfer DM. Functional annotation of rare structural variation in the human brain. BioRxiv doi: https://doi.org/10.1101/711754

## Data Generation

To download/generate annotation files run `bin/make_annotation_files`.

## Annotating CNV
perl bin/annotate-cnv.pl [input_file] > [output_file]

There is a test CNV file in test/test.cnv

## Applying model

The model can be applied using coefficients


## Calculating regulatory disruption score
