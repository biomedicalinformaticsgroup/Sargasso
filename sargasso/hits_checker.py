from collections import namedtuple

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

class HitsChecker:
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
            f._update_hits_info()

        threshold_data = [self._check_thresholds(i, f) for i, f in enumerate(filterers)]
        assignee = self._assign_hits(threshold_data)

        if assignee == REJECTED:
            for filterer in filterers:
                filterer._add_rejected_hits_to_stats()
        elif assignee == AMBIGUOUS:
            for filterer in filterers:
                filterer._add_ambiguous_hits_to_stats()
        else:
            for i, filterer in enumerate(filterers):
                if i == assignee:
                    self.check_and_write_hits_for_read(filterer)
                else:
                    filterer._add_rejected_hits_to_stats()

        for filterer in filterers:
            filterer._clear_hits()

    def check_and_write_hits_for_read(self, filterer):
        if filterer.hits_info is None:
            filterer._update_hits_info()

        if self.check_hits(filterer.hits_info):
            filterer._add_accepted_hits_to_stats()
            filterer._write_hits()
        else:
            filterer._add_rejected_hits_to_stats()

        filterer._clear_hits()

    def check_and_write_hits_for_remaining_reads(self, filterer):
        try:
            while True:
                if filterer.hits_for_read is None:
                    filterer._get_next_read_hits()
                self.check_and_write_hits_for_read(filterer)
        except StopIteration:
            pass

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        return hits_info.get_multimaps() <= self.multimap_thresh and \
            hits_info.get_primary_mismatches() <= \
            round(self.mismatch_thresh * hits_info.get_total_length()) and \
            self._check_cigars(hits_info) != CIGAR_FAIL

    def _assign_hits_standard(self, threshold_data):
        threshold_data = [t for t in threshold_data if not t.violated]

        num_filterers = len(threshold_data)

        if num_filterers == 0:
            return REJECTED
        elif num_filterers == 1:
            return threshold_data[0].index

        min_mismatches = min([m.mismatches for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.mismatches == min_mismatches]

        if len(threshold_data) == 1:
            return threshold_data[0].index

        min_cigar_check = min([m.cigar_check for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.cigar_check == min_cigar_check]

        if len(threshold_data) == 1:
            return threshold_data[0].index

        min_multimaps = min([m.multimaps for m in threshold_data])
        threshold_data = [t for t in threshold_data
                          if t.multimaps == min_multimaps]

        if len(threshold_data) == 1:
            return threshold_data[0].index

        return AMBIGUOUS

    def _assign_hits_reject_multimaps(self, threshold_data):
        if len([t for t in threshold_data if t.multimaps > 1]) > 0:
            return REJECTED

        return self._assign_hits_standard(threshold_data)

    def _check_thresholds(self, index, filterer):
        hits_info = filterer.hits_info
        violated = False

        multimaps = hits_info.get_multimaps()
        if multimaps > self.multimap_thresh:
            violated = True

        mismatches = hits_info.get_primary_mismatches()
        if mismatches > round(self.mismatch_thresh *
                              hits_info.get_total_length()):
            violated = True

        cigar_check = self._check_cigars(hits_info)
        if cigar_check == CIGAR_FAIL:
            violated = True

        return ThresholdData(
            index, violated, multimaps, mismatches, cigar_check)

    def _check_cigars(self, hits_info):
        total_length = hits_info.get_total_length()
        min_match = total_length - round(self.minmatch_thresh * total_length)

        cigars = hits_info.get_primary_cigars()
        response = CIGAR_GOOD

        num_matches = 0;
        for cigar in cigars:
            for operation, length in cigar:
                if operation == CIGAR_OP_MATCH:
                    num_matches += length
                elif operation == CIGAR_OP_REF_INSERTION or \
                        operation == CIGAR_OP_REF_DELETION:
                    response = CIGAR_LESS_GOOD

        if num_matches < min_match:
            return CIGAR_FAIL
        elif num_matches < total_length:
            return CIGAR_LESS_GOOD

        return response
