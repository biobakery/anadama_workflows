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
    def _filter_samples_for_file(sample_group, file_, 
                                 key=lambda val: val[0]):
        return _filter(lambda sample: key(sample) in file_,
                       sample_group)



from .vis import VisualizationPipeline
from .wgs import WGSPipeline
from .sixteen import SixteenSPipeline
