from collections import namedtuple

from sargasso.utlis.factory import Manager


class HitsChecker:
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
        ['index', 'violated', 'multimaps',
         'mismatches', 'cigar_check'])

    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh,
                 reject_multimaps, logger):
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

    def compare_and_write_hits(self, filterers):
        # Compare the hits for a particular read in each species and decide whether
        # the read can be assigned to one species or another, or if it must be
        # rejected as ambiguous
        for f in filterers:
            f.update_hits_info()

        threshold_data = [self._check_thresholds(i, f) for i, f
                          in enumerate(filterers)]

        # # todo remove debug
        # for t in threshold_data:
        #     print(t)

        assignee = self._assign_hits(threshold_data)

        if assignee == self.REJECTED:
            for filterer in filterers:
                filterer.add_rejected_hits_to_stats()
        elif assignee == self.AMBIGUOUS:
            for filterer in filterers:
                filterer.add_ambiguous_hits_to_stats()
        else:
            for i, filterer in enumerate(filterers):
                if i == assignee:
                    self.check_and_write_hits_for_read(filterer)
                else:
                    filterer.add_rejected_hits_to_stats()

        for filterer in filterers:
            filterer.clear_hits()

    def check_and_write_hits_for_read(self, filterer):
        if filterer.hits_info is None:
            filterer.update_hits_info()

        if self.check_hits(filterer.hits_info):
            filterer.add_accepted_hits_to_stats()
            filterer.write_hits()
        else:
            filterer.add_rejected_hits_to_stats()

        filterer.clear_hits()

    def check_and_write_hits_for_remaining_reads(self, filterer):
        try:
            while True:
                if filterer.hits_for_read is None:
                    filterer.get_next_read_hits()
                # # todo remove debug
                #                 # print("Read:{}!!!".format(filterer.hits_for_read[0].qname))
                #                 # print('assigned due to only one competing filterer!')
                self.check_and_write_hits_for_read(filterer)
        except StopIteration:
            pass

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # # todo remove debug multimap
        return hits_info.get_multimaps() <= self.multimap_thresh and \
               hits_info.get_primary_mismatches() <= \
               round(self.mismatch_thresh * hits_info.get_total_length()) and \
               self._check_cigars(hits_info) != self.CIGAR_FAIL

    def _assign_hits_standard(self, threshold_data):
        threshold_data = [t for t in threshold_data if not t.violated]

        num_filterers = len(threshold_data)

        if num_filterers == 0:
            return self.REJECTED
        elif num_filterers == 1:
            # # todo remove debug
            # print('assigned due to only one filter exist!')
            return threshold_data[0].index

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            # # todo remove debug
            # print('assigned due to primary hit min_mismatches!')
            return threshold_data[0].index

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            # # todo remove debug
            # print('assigned due to primart hit CIGAR!')
            return threshold_data[0].index

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            # # todo remove debug multimap
            # print('assigned due to number of multimap!')
            return threshold_data[0].index

        # # todo remove debug
        # print('assigned Ambigous!')
        return self.AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            return self.REJECTED

        return self._assign_hits_standard(threshold_data)

    def _check_thresholds(self, index, filterer):
        hits_info = filterer.hits_info
        violated = False

        multimaps = hits_info.get_multimaps()
        if multimaps > self.multimap_thresh:
            # # todo remove debug multimap
            # print('violated due to multimap!')
            violated = True

        mismatches = hits_info.get_primary_mismatches()
        if mismatches > round(self.mismatch_thresh *
                              hits_info.get_total_length()):
            # # todo remove debug
            # print('violated due to primary mismatches!')
            violated = True

        cigar_check = self._check_cigars(hits_info)
        if cigar_check == self.CIGAR_FAIL:
            # # todo remove debug
            # print('violated due to primary CIGAR!')
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


class RnaseqHitsChecker(HitsChecker):
    pass


class ChipseqHitChecker(HitsChecker):
    pass


class HitsCheckerManager(Manager):
    HITSCHECKERS = {'rnaseq': RnaseqHitsChecker,
                    'chipseq': ChipseqHitChecker}

    @classmethod
    def get(cls, data_type, mismatch_thresh, minmatch_thresh, multimap_thresh,
            reject_multimaps, logger):
        return cls.HITSCHECKERS[data_type](mismatch_thresh, minmatch_thresh, multimap_thresh,
                                           reject_multimaps, logger)

    @classmethod
    def _create(cls, data_type):
        pass
