"""
All anadama_workflows pipelines
"""

import re
import os
from itertools import dropwhile, chain
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


    def _unpack_metadata(self, default=None):
        samples = []
        if type(self.sample_metadata) is str:
            self.sample_metadata = [self.sample_metadata]
        if type(self.sample_metadata[0]) is str:
            with open(self.sample_metadata[0]) as metadata_f:
                samples = list( util.deserialize_map_file(metadata_f) )
        if str(self.sample_metadata[0]).startswith("Sample"):
            samples = self.sample_metadata
                
        if not samples:
            if default:
                self.sample_metadata = default()
            else:
                raise ValueError("Unable to read map.txt file. Empty file.")
        else:
            self.sample_metadata = samples


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


def infer_pairs(list_fnames):
    one_files, two_files, notpairs = _regex_filter(list_fnames)
    if len(one_files) != len(two_files):
        pairs = list()
        notpairs.extend(one_files)
        notpairs.extend(two_files)
    else:
        pairs = zip( sorted(one_files), sorted(two_files) )
    
    return pairs, notpairs


def _regex_filter(list_fnames):
    """Go through each name; group into R1, R2, or singleton based on a
    regular expression that searches for R1 or r2 or R2 or r1.

    """
    regex = re.compile(r'[-._ ][rR]?([12])[-._ ]')
    one, two, notpairs = list(), list(), list()

    matches = zip( list_fnames, map(regex.search, list_fnames))
    for fname, regex_result in matches:
        if not regex_result:
            notpairs.append(fname)
        else:
            if regex_result.group(1) == "1":
                one.append(fname)
            elif regex_result.group(1) == "2":
                two.append(fname)
            else:
                notpairs.append(fname)    

    return one, two, notpairs



def maybe_decompress(raw_seq_files, products_dir):
    if not raw_seq_files:
        idxs, compressed_files = list(), list()
    elif isinstance(raw_seq_files[0], tuple):
        idxs = list(util.which_compressed_idxs(raw_seq_files))
        compressed_files = util.take(raw_seq_files, idxs)
    else:
        packed = zip(*util.filter_compressed(raw_seq_files))
        if packed:
            idxs, compressed_files = packed
        else:
            idxs, compressed_files = list(), list()
    
    tasks = []

    def _decompress(f):
        unz_f = os.path.splitext(f)[0]
        unz_f = util.new_file( os.path.basename(unz_f), basedir=products_dir )
        tasks.append(general.extract(f, unz_f))
        return unz_f

    for idx in idxs:
        if isinstance(idx, tuple):
            raw_seq_files[idx[0]] = list(raw_seq_files[idx[0]])
            compressed = raw_seq_files[idx[0]][idx[1]]
            raw_seq_files[idx[0]][idx[1]] = _decompress(compressed)
        else:
            compressed = raw_seq_files[idx]
            raw_seq_files[idx] = _decompress(compressed)

    return raw_seq_files, tasks


def _to_merged(fname_str, tag="merged", strip_ext=True):
    """tag is put into filename for describing what happened to the
    sequences: `1`merged`` for stitching ``cat`` for concatenation"""
    fname_str = re.sub(r'(.*[-._ ])[rR]?[12]([-._ ].*)', 
                       r'\1%s\2'%(tag), fname_str)
    
    if strip_ext and util.is_compressed(fname_str):
        fname_str = rmext(fname_str)

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
        if guess != "fastq" or util.is_compressed(f):
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
from .rna import RNAPipeline
