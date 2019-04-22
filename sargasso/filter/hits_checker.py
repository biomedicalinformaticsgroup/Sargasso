from collections import namedtuple


class HitsChecker(object):
    REJECTED = -1
    AMBIGUOUS = -2

    CIGAR_GOOD = 0
    CIGAR_LESS_GOOD = 1
    CIGAR_FAIL = 2

    CIGAR_OP_MATCH = 0  # From pysam
    CIGAR_OP_REF_INSERTION = 1  # From pysam
    CIGAR_OP_REF_DELETION = 2  # From pysam
    CIGAR_OP_REF_SKIP = 3  # From pysam

    ThresholdData = namedtuple(
        'ThresholdData',
        ['index', 'violated', 'multimaps', 'mismatches', 'cigar_check'])

    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh,
                 reject_multimaps, logger):
        self.logger = logger
        self.mismatch_thresh = mismatch_thresh / 100.0
        self.minmatch_thresh = minmatch_thresh / 100.0
        self.multimap_thresh = multimap_thresh
        self._assign_hits = self._assign_hits_reject_multimaps \
            if reject_multimaps else self._assign_hits_standard

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

        threshold_data = [self._check_thresholds(i, m) for i, m
                          in enumerate(hits_managers)]

        if __debug__:
            for t in threshold_data:
                self.logger.debug(t)

        assignee = self._assign_hits(threshold_data)

        if assignee == self.REJECTED:
            for hits_manager in hits_managers:
                hits_manager.add_rejected_hits_to_stats()
        elif assignee == self.AMBIGUOUS:
            for hits_manager in hits_managers:
                hits_manager.add_ambiguous_hits_to_stats()
        else:
            for i, hits_manager in enumerate(hits_managers):
                if i == assignee:
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

        if self.check_hits(hits_manager.hits_info):
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

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # @xintodo can we return after each IF to skip the result tests

        violated = False

        if hits_info.get_multimaps() > self.multimap_thresh:
            violated = True
            if __debug__:
                self.logger.debug(
                    '    violated multimap.')

        if hits_info.get_primary_mismatches() > \
                round(self.mismatch_thresh * hits_info.get_total_length()):
            violated = True
            if __debug__:
                self.logger.debug(
                    '    violated primary mismatches.')

        if self._check_cigars(hits_info) == self.CIGAR_FAIL:
            violated = True
            if __debug__:
                self.logger.debug(
                    '    violated primary CIGAR.')

        return not violated

    def _assign_hits_standard(self, threshold_data):
        threshold_data = [t for t in threshold_data if not t.violated]

        num_hits_managers = len(threshold_data)

        if num_hits_managers == 0:
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__:
                self.logger.debug('assigned due to only one filter left after checking threshold!')
            return threshold_data[0].index

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            if __debug__:
                self.logger.debug('assigne due to primary hit min_mismatches!')
            return threshold_data[0].index

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            if __debug__:
                self.logger.debug('assigne due to primart hit CIGAR!')
            return threshold_data[0].index

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            # # todo remove debug multimap
            if __debug__:
                self.logger.debug('assigned due to number of multimap!')
            return threshold_data[0].index

        if __debug__:
            self.logger.debug('assigned due to Ambigous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)

    def _check_thresholds(self, index, hits_manager):

        hits_info = hits_manager.hits_info
        violated = False

        multimaps = hits_info.get_multimaps()
        if multimaps > self.multimap_thresh:
            # # todo remove debug multimap
            if __debug__:
                self.logger.debug('violated due to multimap!')
            violated = True

        mismatches = hits_info.get_primary_mismatches()
        if mismatches > round(self.mismatch_thresh *
                              hits_info.get_total_length()):
            if __debug__:
                self.logger.debug('violated due to primary mismatches!')
            violated = True

        cigar_check = self._check_cigars(hits_info)
        if cigar_check == self.CIGAR_FAIL:
            if __debug__:
                self.logger.debug('violated due to primary CIGAR!')
            violated = True

        return self.ThresholdData(
            index, violated, multimaps, mismatches, cigar_check)

    def _check_cigars(self, hits_info):
        total_length = hits_info.get_total_length()
        min_match = total_length - round(self.minmatch_thresh * total_length)

        cigars = hits_info.get_primary_cigars()
        response = self.CIGAR_GOOD

        num_matches = 0
        for cigar in cigars:
            for operation, length in cigar:
                if operation == self.CIGAR_OP_MATCH:
                    num_matches += length
                elif operation == self.CIGAR_OP_REF_INSERTION or \
                        operation == self.CIGAR_OP_REF_DELETION:
                    response = self.CIGAR_LESS_GOOD

        if num_matches < min_match:
            return self.CIGAR_FAIL
        elif num_matches < total_length:
            return self.CIGAR_LESS_GOOD

        return response


class RnaSeqHitsChecker(HitsChecker):
    pass


class DnaSeqHitsChecker(HitsChecker):
    pass


class BisulfiteHitsChecker(HitsChecker):

    ThresholdData = namedtuple(
        'ThresholdData',
        ['index', 'violated', 'multimaps', 'mismatches', 'cigar_check', 'is_ambig'])

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        violated = not super(self.__class__, self).check_hits(hits_info)

        if hits_info.get_is_ambig_hit():
            if __debug__:
                self.logger.debug('    violated ambig.')
            violated = True

        return not violated

    def _check_thresholds(self, index, hits_manager):

        hits_info = hits_manager.hits_info

        if __debug__:
            self.logger.debug("  species {}:".format(hits_manager.species_id))

        violated = not self.check_hits(hits_info)

        return self.ThresholdData(
            index, violated, hits_info.get_multimaps(), hits_info.get_primary_mismatches(),
            hits_info.get_primary_cigars(), hits_info.get_is_ambig_hit())

    def _assign_hits_standard(self, threshold_data):

        unique = [not t.is_ambig for t in threshold_data]
        # if any(unique):
        #     ## we need to keep all the threshold data, is any of it has is_ambig=False
        #     threshold_data=threshold_data
        #     if any([t.is_ambig for t in threshold_data]):
        #         print('debug break point')
        # else:
        #     ## all hits are ambig, skip
        #     threshold_data=[]

        if all(unique):
            threshold_data = [t for t in threshold_data if not t.violated]
        else:
            ## all hits are ambig, skip
            threshold_data=[]

        num_hits_managers = len(threshold_data)

        if num_hits_managers == 0:
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__:
                self.logger.debug('assigned due to only one filter left after checking threshold!')
            return threshold_data[0].index

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            if __debug__:
                self.logger.debug('assigne due to primary hit min_mismatches!')
            return threshold_data[0].index

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            if __debug__:
                self.logger.debug('assigne due to primart hit CIGAR!')
            return threshold_data[0].index

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            # # todo remove debug multimap
            if __debug__:
                self.logger.debug('assigned due to number of multimap!')
            return threshold_data[0].index

        if __debug__:
            self.logger.debug('assigned due to Ambigous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)
