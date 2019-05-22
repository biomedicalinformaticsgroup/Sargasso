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
        self.threshold_data_cls = ThresholdData

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
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__: self.logger.debug('assigned due to only one filter left after checking threshold!')
            return threshold_data[0].species_id

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            if __debug__: self.logger.debug('assigne due to primary hit min_mismatches!')
            return threshold_data[0].species_id

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            if __debug__: self.logger.debug('assigne due to primart hit CIGAR!')
            return threshold_data[0].species_id

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            # # todo remove debug multimap
            if __debug__:
                self.logger.debug('assigned due to number of multimap!')
            return threshold_data[0].species_id

        if __debug__:
            self.logger.debug('assigned due to Ambigous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)

    def _check_thresholds(self, hits_manager):
        return self.threshold_data_cls(hits_manager, self.mismatch_thresh, self.minmatch_thresh, self.multimap_thresh,
                                       self.logger)


class RnaSeqHitsChecker(HitsChecker):
    pass


class DnaSeqHitsChecker(HitsChecker):
    pass


class BisulfiteHitsChecker(HitsChecker):

    def _check_thresholds(self, hits_manager):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # return True if valid
        threshold_data = super(self.__class__, self)._check_thresholds(hits_manager)

        violated_ambig = False

        is_ambig = hits_manager.hits_info.get_is_ambig_hit()
        if is_ambig:
            if __debug__: self.logger.debug('    violated ambig.')
            violated_ambig = True

        threshold_data.violated = threshold_data.violated or violated_ambig
        threshold_data.is_ambig = is_ambig
        threshold_data.violated_ambig = violated_ambig

        return threshold_data

    def _assign_hits_standard(self, threshold_data):
        ## we need to keep all the threshold data, if any of it has is_ambig=False
        ## in order to apply the following logic for the separation
        ## check here for the logic https://github.com/statbio/Sargasso/wiki/Bisulfite-sequencing
        if all([t.is_ambig for t in threshold_data]):
            ## all hits are ambig, skip
            threshold_data = []

        num_hits_managers = len(threshold_data)

        if num_hits_managers == 0:
            return self.REJECTED
        elif num_hits_managers == 1:
            if __debug__: self.logger.debug('assigned due to only one filter left after checking threshold!')
            return threshold_data[0].index

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            ## if the better read is from the ambig species, we reject
            if threshold_data[0].is_ambig:
                return self.REJECTED
            else:
                return threshold_data[0].index

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            ## if the better read is from the ambig species, we reject
            if threshold_data[0].is_ambig:
                return self.REJECTED
            else:
                return threshold_data[0].index

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
        #         return threshold_data[0].index
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        ## for bisulfite, the reject_multimaps means the reads is ambigours
        if any([t for t in threshold_data if t.is_ambig]):
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)


class ThresholdData(object):
    CIGAR_GOOD = 0
    CIGAR_LESS_GOOD = 1
    CIGAR_FAIL = 2

    CIGAR_OP_MATCH = 0  # From pysam
    CIGAR_OP_REF_INSERTION = 1  # From pysam
    CIGAR_OP_REF_DELETION = 2  # From pysam
    CIGAR_OP_REF_SKIP = 3  # From pysam

    def __init__(self, hits_manager, mismatch_thresh, minmatch_thresh, multimap_thresh, logger):

        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # return True if valid

        self.species_id = hits_manager.species_id

        hits_info = hits_manager.hits_info

        self.violated_multimaps = False
        self.violated_mismatches = False
        self.violated_cigar_check = False

        self.multimaps = hits_info.get_multimaps()
        if self.multimaps > multimap_thresh:
            if __debug__: logger.debug('violated due to multimap!')
            self.violated_multimaps = True

        self.mismatches = hits_info.get_primary_mismatches()
        if self.mismatches > round(mismatch_thresh *
                                   hits_info.get_total_length()):
            if __debug__: logger.debug('violated due to primary mismatches!')
            self.violated_mismatches = True

        self.cigar_check = self._check_cigars(hits_info, minmatch_thresh)
        if self.cigar_check == self.CIGAR_FAIL:
            if __debug__: logger.debug('violated due to primary CIGAR!')
            self.violated_cigar_check = True

        self.violated = any([self.violated_multimaps, self.violated_mismatches, self.violated_cigar_check])

    def _check_cigars(self, hits_info, minmatch_thresh):
        total_length = hits_info.get_total_length()
        min_match = total_length - round(minmatch_thresh * total_length)

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

# class BisulfiteThresholdData(object):
# #     def __init__(self, index, violated, multimaps, mismatches, cigar_check,
# #                  violated_multimaps, violated_mismatches, violated_cigar_check):
# #         self.index = index
# #         self.violated = violated
# #         self.multimaps = multimaps
# #         self.mismatches = mismatches
# #         self.cigar_check = cigar_check
# #         self.violated_multimaps = violated_multimaps
# #         self.violated_mismatches = violated_mismatches
# #         self.violated_cigar_check = violated_cigar_check
# #         self.is_ambig = is_ambig
# #         self.violated_ambig = violated_ambig
