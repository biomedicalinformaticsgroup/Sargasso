from . import hits_checker
from . import hits_info
from . import samutils
from . import separation_stats


class Filterer(object):
    def __init__(self, species_id, input_bam, output_bam, h_check, logger):
        self.species_id = species_id
        self.stats = separation_stats.SeparationStats(species_id)

        input_hits = samutils.open_samfile_for_read(input_bam)
        self.output_bam = samutils.open_samfile_for_write(
            output_bam, input_hits)

        self.hits_generator = samutils.hits_generator(input_hits)
        self.hits_for_read = None
        self.hits_info = None
        self.count = 0
        self.hits_checker = h_check
        self.logger = logger

    def get_next_read_name(self):
        if self.hits_for_read is None:
            self._get_next_read_hits()
            self.count += 1
            if self.count % 1000000 == 0:
                self.logger.debug("Read {n} reads from species {s}".format(
                    n=self.count, s=self.species_id))

        return self.hits_for_read[0].query_name

    def compare_and_write_hits(self, other_filterer):
        # Compare the hits for a particular read in each species and decide whether
        # the read can be assigned to one species or another, or if it must be
        # rejected as ambiguous
        self._update_hits_info()
        other_filterer._update_hits_info()
        # check read length to decide a weighting for the filter criteria
        self.hits_checker._calculate_weighting(self.hits_info)
        other_filterer.hits_checker._calculate_weighting(self.hits_info)

        assignment = self.hits_checker.assign_hits(self.hits_info, other_filterer.hits_info)

        if assignment == hits_checker.ASSIGNED_TO_SPECIES_ONE:
            other_filterer._add_rejected_hits_to_stats()
            self.check_and_write_hits_for_read()
        elif assignment == hits_checker.ASSIGNED_TO_SPECIES_TWO:
            self._add_rejected_hits_to_stats()
            other_filterer.check_and_write_hits_for_read()
        elif assignment == hits_checker.ASSIGNED_TO_NEITHER_REJECTED:
            self._add_rejected_hits_to_stats()
            other_filterer._add_rejected_hits_to_stats()
        else:
            self._add_ambiguous_hits_to_stats()
            other_filterer._add_ambiguous_hits_to_stats()

        self._clear_hits()
        other_filterer._clear_hits()

    def check_and_write_hits_for_read(self):
        if self.hits_info is None:
            self._update_hits_info()

        if self.hits_checker.check_hits(self.hits_info):
            self._add_accepted_hits_to_stats()
            self._write_hits()
        else:
            self._add_rejected_hits_to_stats()

        self._clear_hits()

    def check_and_write_hits_for_remaining_reads(self):
        try:
            while True:
                # TODO: is this potentially skipping a read?
                self._get_next_read_hits()
                self.check_and_write_hits_for_read()
        except StopIteration:
            pass

    def log_stats(self):
        self.logger.info(self.stats)

    def _get_next_read_hits(self):
        self.hits_for_read = self.hits_generator.next()
        self.hits_info = None

    def _update_hits_info(self):
        self.hits_info = hits_info.HitsInfo(self.hits_for_read)

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
