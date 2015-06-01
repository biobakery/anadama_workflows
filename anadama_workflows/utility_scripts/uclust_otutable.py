"""Goes from a usearch mapping results file to an OTU table, 
-> samples \/ OTUs 

Assumes that search query sequences are from qiime-formatted
demultiplexed sequence files. That is, the sequence label is composed
of the sampleid along with other items, separated by underscores '_'.
The first item before the underscore is assumed to be the sample id.

"""

import re
import sys
from collections import defaultdict

def parsetarget(target_str):
    return re.search(r'OTU_(\d+)', target_str).group(1)

def fields(uc_fname):
    with open(uc_fname) as f:
        for line in f:
            if line.startswith("H"):
                fields = line.split('\t')
                query, target = fields[8], fields[9]
                sample_id = query.split("_", 1)[0]
                otu_id = parsetarget(target)
                yield sample_id, otu_id

def get(d, items, default=0):
    return [d.get(item, default) for item in items]

def output(d, sample_ids):
    sample_ids = list(sample_ids)
    print "\t".join(["OTUId"]+list(sample_ids))
    for otu_id, hits_dict in d.iteritems():
        hits = get(hits_dict, sample_ids, default=0)
        print "\t".join([otu_id]+map(str, hits))

    
def main(uc_fname=None):
    if not uc_fname:
        uc_fname = sys.argv[1]
    table_dict = defaultdict(lambda:defaultdict(int))
    sample_ids = set()
    for sample_id, otu_id in fields(uc_fname):
        sample_ids.add(sample_id)
        table_dict[otu_id][sample_id] += 1

    output(table_dict, sample_ids)

if __name__ == "__main__":
    ret = main(*sys.argv[1:])
    sys.exit(ret)
