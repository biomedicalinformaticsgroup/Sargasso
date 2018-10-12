from factory import Manager
from samutils import SamUtilsManager
from . import hits_info


class Filterer(object):
    def __init__(self, data_type, species_id, input_bam, output_bam, logger):
        self.data_type = data_type
        self.species_id = species_id
        self.stats = SeparationStats(species_id)
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
            self.get_next_read_hits()
            self.count += 1
            if self.count % 1000000 == 0:
                self.logger.debug("Read {n} reads from species {s}".format(
                    n=self.count, s=self.species_id))

        return self.hits_for_read[0].query_name

    def log_stats(self):
        self.logger.info(self.stats)

    def get_next_read_hits(self):
        self.hits_for_read = self.hits_generator.next()
        self.hits_info = None

    def update_hits_info(self):
        self.hits_info = hits_info.HitsInfoManager.get(self.data_type, self.hits_for_read)

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


class RnaseqFilterer(Filterer):
    pass


class ChipseqFilterer(Filterer):
    pass


class FilterManager(Manager):
    FILTERS = {'rnaseq': RnaseqFilterer,
               'chipseq': ChipseqFilterer}

    @staticmethod
    def get(data_type, species_id, input_bam, output_bam, logger):
        return FilterManager.FILTERS[data_type](data_type, species_id, input_bam, output_bam, logger)


class SeparationStats:
    def __init__(self, species_id):
        self.name = "Species {n}".format(n=species_id)
        self.hits_written = 0
        self.reads_written = 0
        self.hits_rejected = 0
        self.reads_rejected = 0
        self.hits_ambiguous = 0
        self.reads_ambiguous = 0

    def accepted_hits(self, hits):
        self.hits_written += len(hits)
        self.reads_written += 1

    def rejected_hits(self, hits):
        self.hits_rejected += len(hits)
        self.reads_rejected += 1

    def ambiguous_hits(self, hits):
        self.hits_ambiguous += len(hits)
        self.reads_ambiguous += 1

    def __str__(self):
        return ("{n}: wrote {f} filtered hits for {fr} reads; {r} hits for " +
                "{rr} reads were rejected outright, and {a} hits for " +
                "{ar} reads were rejected as ambiguous.").format(
            n=self.name,
            f=self.hits_written, fr=self.reads_written,
            r=self.hits_rejected, rr=self.reads_rejected,
            a=self.hits_ambiguous, ar=self.reads_ambiguous)
