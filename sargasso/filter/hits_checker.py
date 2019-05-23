from sargasso.filter.threshold_checker import *


class HitsChecker(object):
    REJECTED = -1
    AMBIGUOUS = -2

    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh,
                 reject_multimaps, logger):
        self.logger = logger
        self.mismatch_thresh = mismatch_thresh / 100.0
        self.minmatch_thresh = minmatch_thresh / 100.0
        self.multimap_thresh = multimap_thresh
        self._assign_hits = self._assign_hits_reject_multimaps \
            if reject_multimaps else self._assign_hits_standard
        self.threshold_checker_cls = ThresholdChecker

        logger.debug(("PARAMS: mismatch - {mism}, minmatch - {minm}, " +
                      "multimap - {mult}").format(
            mism=self.mismatch_thresh,
            minm=self.minmatch_thresh,
            mult=self.multimap_thresh))

    def compare_and_write_hits(self, hits_managers):
        # Compare the hits for a particular read in each species and decide whether
        # the read can be assigned to one species or another, or if it must be
        # rejected as ambiguous
        for m in hits_managers:
            m.update_hits_info()

        threshold_data = [self._check_thresholds(m) for i, m in enumerate(hits_managers)]

        assignee = self._assign_hits(threshold_data)

        if assignee == self.REJECTED:
            for hits_manager in hits_managers:
                hits_manager.add_rejected_hits_to_stats()
        elif assignee == self.AMBIGUOUS:
            for hits_manager in hits_managers:
                hits_manager.add_ambiguous_hits_to_stats()
        else:
            for i, hits_manager in enumerate(hits_managers):
                if hits_manager.species_id == assignee:
                    hits_manager.add_accepted_hits_to_stats()
                    hits_manager.write_hits()
                    # self.check_and_write_hits_for_read(hits_manager)
                else:
                    hits_manager.add_rejected_hits_to_stats()

        for hits_manager in hits_managers:
            hits_manager.clear_hits()

    def check_and_write_hits_for_read(self, hits_manager):
        if hits_manager.hits_info is None:
            hits_manager.update_hits_info()

        if self.check_hits(hits_manager):
            hits_manager.add_accepted_hits_to_stats()
            hits_manager.write_hits()
        else:
            hits_manager.add_rejected_hits_to_stats()

        hits_manager.clear_hits()

    def check_and_write_hits_for_remaining_reads(self, hits_manager):
        try:
            while True:
                if hits_manager.hits_for_read is None:
                    hits_manager.get_next_read_hits()
                self.check_and_write_hits_for_read(hits_manager)
        except StopIteration:
            pass

    def check_hits(self, hits_manager):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # return True if valid
        threshold_data = self._check_thresholds(hits_manager)
        return not threshold_data.violated

    def _assign_hits_standard(self, threshold_data):

        threshold_data = [t for t in threshold_data if not t.violated]

        num_hits_managers = len(threshold_data)

        if num_hits_managers == 0:
            if __debug__: self.logger.debug('Reject: All violated!')
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__: self.logger.debug('Assign {}: only valid threshold'.format(threshold_data[0].species_id))
            return threshold_data[0].species_id

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            if __debug__: self.logger.debug('Assign {}: min_mismatches!'.format(threshold_data[0].species_id))
            return threshold_data[0].species_id

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            if __debug__: self.logger.debug('Assign {}: CIGAR!'.format(threshold_data[0].species_id))
            return threshold_data[0].species_id

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            if __debug__: self.logger.debug('Assign {}: number of multimap!'.format(threshold_data[0].species_id))
            return threshold_data[0].species_id

        if __debug__:
            self.logger.debug('Reject: Ambigous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            if __debug__: self.logger.debug('Reject: is multimap!')
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)

    def _check_thresholds(self, hits_manager):
        thc = self.threshold_checker_cls(hits_manager,
                                         self.mismatch_thresh,
                                         self.minmatch_thresh,
                                         self.multimap_thresh)
        if __debug__:
            self.logger.debug(thc.to_debug_string())
        return thc


class RnaSeqHitsChecker(HitsChecker):
    pass


class DnaSeqHitsChecker(HitsChecker):
    pass


class BisulfiteHitsChecker(HitsChecker):
    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh,
                 reject_multimaps, logger):
        super(self.__class__, self).__init__(mismatch_thresh, minmatch_thresh, multimap_thresh,
                                             reject_multimaps, logger)
        self.threshold_checker_cls = BisulfiteThresholdChecker

    def _assign_hits_standard(self, threshold_data):
        ## we need to keep all the threshold data, if any of it has is_ambig=False
        ## in order to apply the following logic for the separation
        ## check here for the logic https://github.com/statbio/Sargasso/wiki/Bisulfite-sequencing
        if all([t.is_ambig for t in threshold_data]):
            ## all hits are ambig, skip
            threshold_data = []

        num_hits_managers = len(threshold_data)

        if num_hits_managers == 0:
            if __debug__: self.logger.debug('Reject: All ambig!')
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__: self.logger.debug('Assign {}: only valid threshold'.format(threshold_data[0].species_id))
            return threshold_data[0].species_id

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            ## if the better read is from the ambig species, we reject
            if threshold_data[0].is_ambig:
                if __debug__: self.logger.debug(
                    'Reject: ambig has a better min_mismatches!'.format(threshold_data[0].species_id))
                return self.REJECTED
            else:
                if __debug__: self.logger.debug('Assign {}: min_mismatches!'.format(threshold_data[0].species_id))
                return threshold_data[0].species_id

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            ## if the better read is from the ambig species, we reject
            if threshold_data[0].is_ambig:
                if __debug__: self.logger.debug(
                    'Reject: ambig has a better CIGAR!'.format(threshold_data[0].species_id))
                return self.REJECTED
            else:
                if __debug__: self.logger.debug('Assign {}: CIGAR!'.format(threshold_data[0].species_id))
                return threshold_data[0].species_id

        ## bismark only return A best hit, so the number of multimap will always be 1.
        ## The actural number of multimap is unknown.
        ## Thus we do not need to check multimap for bismark output
        # min_multimaps = min([m.multimaps for m in threshold_data])
        # threshold_data = [t for t in threshold_data
        #                   if t.multimaps == min_multimaps]
        #
        # if len(threshold_data) == 1:
        #     ## if the better read is from the ambig species, we reject
        #     if threshold_data[0].is_ambig:
        #         return self.REJECTED
        #     else:
        #         if __debug__:
        #             self.logger.debug('assigned due to number of multimap!')
        #         return threshold_data[0].species_id
        if __debug__: self.logger.debug('Ambiguous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        ## for bisulfite, the reject_multimaps means the reads is ambigours
        ## This is check first by the _assign_hits_standard, thus no need to check here
        return self._assign_hits_standard(threshold_data)


