import sargasso.utils.samutils as su

from sargasso.filter import hits_info
from sargasso.filter.separation_stats import SeparationStats


class HitsManager(object):
    def __init__(
        self, hits_info_cls, species_id, input_bam, output_bam, logger):

        self.hits_info_cls = hits_info_cls
        self.species_id = species_id
        self.stats = SeparationStats(species_id)

        input_hits = su.open_samfile_for_read(input_bam)
        self.output_bam = su.open_samfile_for_write(output_bam, input_hits)

        self.hits_generator = su.hits_generator(input_hits)
        self.hits_for_read = None
        self.hits_info = None
        self.count = 0
        self.logger = logger

    def get_next_read_name(self):
        if self.hits_for_read is None:
            self.get_next_read_hits()
            self.count += 1
            if self.count % 1000000 == 0:
                self.logger.debug("Read {n} reads from species {s}".format(
                    n=self.count, s=self.species_id))

        return self.hits_for_read[0].query_name

    def log_stats(self):
        self.logger.info(self.stats)

    def get_next_read_hits(self):
        self.hits_for_read = next(self.hits_generator)
        self.hits_info = None

    def update_hits_info(self):
        self.hits_info = self.hits_info_cls(self.hits_for_read)

    def write_hits(self):
        for hit in self.hits_for_read:
            self.output_bam.write(hit)

    def clear_hits(self):
        self.hits_for_read = None
        self.hits_info = None

    def add_accepted_hits_to_stats(self):
        self.stats.accepted_hits(self.hits_for_read)

    def add_rejected_hits_to_stats(self):
        self.stats.rejected_hits(self.hits_for_read)

    def add_ambiguous_hits_to_stats(self):
        self.stats.ambiguous_hits(self.hits_for_read)


class RnaSeqHitsManager(HitsManager):
    def __init__(self, species_id, input_bam, output_bam, logger):
        HitsManager.__init__(
            self, hits_info.RnaSeqHitsInfo, species_id,
            input_bam, output_bam, logger)


class DnaSeqHitsManager(HitsManager):
    def __init__(self, species_id, input_bam, output_bam, logger):
        HitsManager.__init__(
            self, hits_info.DnaSeqHitsInfo, species_id,
            input_bam, output_bam, logger)
