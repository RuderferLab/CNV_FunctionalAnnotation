###
#   name      | mary lauren benton
#   python_v  | python3
#   created   | 2018.09.06
#   updated   | 2018.11.19
#             | 2018.11.26
#
###

import pandas as pd
import numpy  as np

from scipy import stats
from pybedtools import BedTool


def filter_chipseq_files(infile):
    df = pd.read_table(infile, sep='\t', header=None,
                       names=['chrom', 'start', 'end', 'name', 'score',
                              'strand', 'signal', 'p', 'q', 'peak'])
    df['signal_zscore'] = stats.zscore(df.signal)
    fltr_df = df[df.signal_zscore > 1.64]
    return BedTool.from_dataframe(fltr_df)


# make promoter file | general promoter set
genes = BedTool('rnaseq_ensembl_genemodels_filter.bed')
tss_window = genes.flank(l=2000, r=0, s=True, genome='hg19')
tss_window.saveas('hg19_tss_2kb.bed')

# make promoter file | requires overlap with high signal H3K27ac/H3K4me3
h3k27ac = filter_chipseq_files('E073-H3K27ac.narrowPeak.gz')
h3k4me3 = filter_chipseq_files('E073-H3K4me3.narrowPeak.gz')
tss_high_conf = (tss_window.intersect(h3k27ac, u=True)).intersect(h3k4me3, u=True)
tss_high_conf.saveas('hg19_tss_high_conf.bed')

