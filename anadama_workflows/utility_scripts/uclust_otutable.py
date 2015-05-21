import re
import sys
from collections import defaultdict

def parsetarget(target_str):
    otu_id = re.search(r'OTU_(\d+)', target_str).group(1)
    match = re.search(r'size=(\d+)', target_str)
    size = 1
    if match:
        size = int(match.group(1))
    return otu_id, size

def fields(uc_fname):
    with open(uc_fname) as f:
        for line in f:
            if line.startswith("H"):
                fields = line.split('\t')
                query, target = fields[8], fields[9]
                sample_id = query.split("_", 1)[0]
                otu_id, size = parsetarget(target)
                yield sample_id, otu_id, size

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
    for sample_id, otu_id, size in fields(uc_fname):
        sample_ids.add(sample_id)
        table_dict[otu_id][sample_id] += size

    output(table_dict, sample_ids)

if __name__ == "__main__":
    ret = main(*sys.argv[1:])
    sys.exit(ret)
