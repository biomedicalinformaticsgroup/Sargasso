from . import hits_info
from . import separation_stats
from samutils import SamUtilsManager



class Filterer(object):
    def __init__(self, data_type, species_id, input_bam, output_bam, logger):
        self.data_type = data_type
        self.species_id = species_id
        self.stats = separation_stats.SeparationStats(species_id)
        self.samutils = SamUtilsManager.get(data_type)

        input_hits = self.samutils.open_samfile_for_read(input_bam)
        self.output_bam = self.samutils.open_samfile_for_write(
            output_bam, input_hits)

        self.hits_generator = self.samutils.hits_generator(input_hits)
        self.hits_for_read = None
        self.hits_info = None
        self.count = 0
        self.logger = logger

    def get_next_read_name(self):
        if self.hits_for_read is None:
            self._get_next_read_hits()
            self.count += 1
            if self.count % 1000000 == 0:
                self.logger.debug("Read {n} reads from species {s}".format(
                    n=self.count, s=self.species_id))

        return self.hits_for_read[0].query_name

    def log_stats(self):
        self.logger.info(self.stats)

    def _get_next_read_hits(self):
        self.hits_for_read = next(self.hits_generator)
        self.hits_info = None

    def _update_hits_info(self):
        self.hits_info = hits_info.HitsInfoManager.get(self.data_type, self.hits_for_read)

    def _write_hits(self):
        for hit in self.hits_for_read:
            self.output_bam.write(hit)

    def _clear_hits(self):
        self.hits_for_read = None
        self.hits_info = None

    def _add_accepted_hits_to_stats(self):
        self.stats.accepted_hits(self.hits_for_read)

    def _add_rejected_hits_to_stats(self):
        self.stats.rejected_hits(self.hits_for_read)

    def _add_ambiguous_hits_to_stats(self):
        self.stats.ambiguous_hits(self.hits_for_read)

class RnaseqFilterer(Filterer):
    pass
class ChipseqFilterer(Filterer):
    pass



from factory import Manager
class FilterManager(Manager):
    FILTERS = {'rnaseq':RnaseqFilterer,
               'chipseq': ChipseqFilterer}
    @staticmethod
    def get(data_type, species_id, input_bam, output_bam, logger):
        return FilterManager.FILTERS[data_type](data_type, species_id, input_bam, output_bam, logger)
