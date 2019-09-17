# CNV_FunctionalAnnotation (All currently done in hg19)
Details to annotate CNV as described in Han et al. BioRxiv 2019

Required python3 dependencies: defaultdict, pybedtools, numpy, pandas, scipy

## Data Generation

To download/generate annotation files run `bin/make_annotation_files`.

## Annotating CNV
perl bin/annotate-cnv.pl [input_file] > [output_file]

There is a test CNV file in test/test.cnv


