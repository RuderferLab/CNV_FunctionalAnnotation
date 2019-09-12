###
#   name      | mary lauren benton
#   conda_env | enh_gain-loss
#   created   | 2019.04.04
#
#   this script will process the files downloaded in accessions.txt and merge
#   chip-seq from the same tissues
#   sorted and merged files saved to ./ctcf/encode_tissues/
###


from collections import defaultdict
from pybedtools import BedTool


name_dict = defaultdict(list)

with open('accessions_mapped_filenames.txt') as infile:
    next(infile)
    for line in infile:
        data = line.strip().split('\t')
        tissue_type = data[1].replace(' ', '_').replace('\'', '').lower()
        name_dict[tissue_type].append(data[0])

for tissue in name_dict.keys():
    num_files = len(name_dict[tissue])
    for idx, bed in enumerate(name_dict[tissue]):
        if idx == 0:
            a = BedTool(bed + '.bed')
        else:
            a = a.cat(BedTool(bed + '.bed'), postmerge=True)
    a.sort().merge().saveas(tissue + '.bed')

