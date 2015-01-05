import re
import os
from itertools import dropwhile, chain, izip_longest
from anadama import util

from .. import general


def _filter(func, iterable):
    try:
        filtered = filter(func, iterable)
    except AttributeError:
        return iterable
    return filtered or iterable


class SampleFilterMixin(object):

    @staticmethod
    def _filter_files_for_sample(files_list, sample_group, 
                                 key=lambda val: getattr(val, "Run_accession")):
        return _filter(lambda f: any(key(s) in f for s in sample_group),
                       files_list )

    @staticmethod
    def _filter_pairs_for_sample(files_list, sample_group, 
                                 key=lambda val: getattr(val, "Run_accession")):
        return _filter(lambda f: any(key(s) in f[0] for s in sample_group),
                       files_list )

    @staticmethod
    def _filter_samples_for_file(sample_group, file_, 
                                 key=lambda val: val[0]):
        return _filter(lambda sample: key(sample) in file_,
                       sample_group)


class SampleMetadataMixin(object):

    @property
    def _inferred_sample_metadata_fname(self):
        isnot_str = lambda val: type(val) is not str
        nonstr_products = dropwhile(isnot_str, 
                                    chain(*self.products.itervalues()))
        try:
            base_input = nonstr_products.next()
        except StopIteration:
            raise ValueError("Unable to infer map.txt file location"
                             " because pipeline inputs are empty")

        dir_ = os.path.dirname(base_input)
        return os.path.join(dir_, "map.txt")


    def _get_or_create_sample_metadata(self):
        if type(self.sample_metadata) is not str:
            sample_metadata_fname = self._inferred_sample_metadata_fname
            try:
                util.serialize_map_file(self.sample_metadata, 
                                        sample_metadata_fname)
            except IndexError as e:
                raise ValueError("The provided sample metadata is not in list"
                                 " format, nor is it a string. Sample_metadata"
                                 " should either be a string for a map.txt"
                                 " filename or a list of namedtuples"
                                 " representing the sample metadata")
            return sample_metadata_fname
        else:
            if not os.path.exists(self.sample_metadata):
                raise ValueError("The provided sample metadata file "
                                 "does not exist: "+self.sample_metadata)
            else:
                return self.sample_metadata



def maybe_stitch(maybe_pairs, products_dir, barcode_files=list()):
    pairs, singles = split_pairs(maybe_pairs)
    tasks = list()
    barcodes = list()
    for pair, maybe_barcode in izip_longest(pairs, barcode_files):
        (forward, reverse), maybe_tasks = maybe_convert_to_fastq(
            pair, products_dir)
        tasks.extend(maybe_tasks)
        output = util.new_file( 
            _to_merged(forward),
            basedir=products_dir 
        )
        singles.append(output)
        tasks.append( general.fastq_join(forward, reverse, output) )
        if maybe_barcode:
            filtered_barcode = util.new_file(
                addtag(maybe_barcode, "filtered"),
                basedir=products_dir
            )
            pairtask = general.sequence_pair(
                maybe_barcode, output,
                outfname1=filtered_barcode,
                options={"inner_join": "right"}
            )
            barcodes.append(filtered_barcode)
            tasks.append(pairtask)

    return singles, barcodes, tasks

def maybe_decompress(raw_seq_files):
    if isinstance(raw_seq_files[0], tuple):
        idxs = list(util.which_compressed_idxs(raw_seq_files))
        compressed_files = util.take(raw_seq_files, idxs)
    else:
        idxs, compressed_files = zip(*util.filter_compressed(raw_seq_files))
    
    comp_raw_seq_files = []

    for idx in idxs:
        if isinstance(idx, tuple):
            raw_seq_files[idx[0]] = list(raw_seq_files[idx[0]])
            comp_raw_seq_file = raw_seq_files[idx[0]][idx[1]]
            comp_raw_seq_files.append(comp_raw_seq_file)
            raw_seq_files[idx[0]][idx[1]] = os.path.splitext(
                comp_raw_seq_file)[0]
        else:
            comp_raw_seq_file = raw_seq_files[idx]
            comp_raw_seq_files.append(comp_raw_seq_file)
            raw_seq_files[idx] = os.path.splitext(comp_raw_seq_file)[0]

    return raw_seq_files, comp_raw_seq_files


def _to_merged(fname_str):
    fname_str = re.sub(r'(.*)[rR]1(.*)', r'\1merged\2', fname_str)
    if fname_str.endswith(".gz"):
        fname_str = fname_str[:-3]

    return fname_str


def split_pairs(maybe_pairs):
    pairs, notpairs = list(), list()
    for item in maybe_pairs:
        if (type(item) is list or type(item) is tuple) \
            and len(item) > 1:
            pairs.append(item)
        else:
            notpairs.append(item)

    return pairs, notpairs


def maybe_convert_to_fastq(fnames, products_dir):
    new_fnames, tasks = list(), list()
    for f in fnames:
        guess = util.guess_seq_filetype(f)
        if guess != "fastq" or f.endswith(".bz2"):
            fastq_file = util.new_file(f+".fastq", basedir=products_dir)
            new_fnames.append(fastq_file)
            tasks.append(
                general.sequence_convert([f], fastq_file)
            )
        else:
            new_fnames.append(f)

    return new_fnames, tasks
            

from .vis import VisualizationPipeline
from .wgs import WGSPipeline
from .sixteen import SixteenSPipeline
