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
                 reject_multimaps, overhang_threshold, logger):
        self.mismatch_thresh = mismatch_thresh / 100.0
        self.minmatch_thresh = minmatch_thresh / 100.0
        self.multimap_thresh = multimap_thresh
        self._assign_hits = self._assign_hits_reject_multimaps \
            if reject_multimaps else self._assign_hits_standard
        self.overhang_threshold = overhang_threshold

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
        #total read lenth
        read_length=hits_info.get_total_length()

        min_match = read_length - round(self.minmatch_thresh * read_length)

        cigars = hits_info.get_primary_cigars()

        # when allowing other params, this is no longer a t/f scenario;
        # e.g. when clipping is allowed a cigar without clipping should score
        # better than one with even though both are allowed
        # graded_response = CIGAR_GOOD

        graded_response = self._check_cigar(cigars, min_match, read_length)

        # for cigar in cigar:
        #     graded_response = self._check_cigar(
        #         cigar, graded_response, min_match, read_length)
        #     if graded_response == CIGAR_FAIL:
        #         return graded_response

        return graded_response


    # def _check_cigar(self, cigar, current_response, min_match, length):
    #     graded_response = current_response
    #
    #     num_matches = self._get_cigar_total_match_length(cigar)
    #     if num_matches < min_match:
    #         return CIGAR_FAIL
    #     elif num_matches < length:
    #         graded_response = CIGAR_LESS_GOOD
    #
    #     if self._get_cigar_contains_intron(cigar) and \
    #             not self._check_min_contiguous_match(cigar):
    #         return CIGAR_FAIL
    #
    #     return graded_response

    def _check_cigar(self, cigars, min_match, length):
        graded_response = CIGAR_GOOD

        num_matches = 0;
        for cigar in cigars:
            num_matches += self._get_cigar_total_match_length(cigar)

        if num_matches < min_match:
            return CIGAR_FAIL
        elif num_matches < length:
            graded_response = CIGAR_LESS_GOOD

        for cigar in cigars:
            if self._get_cigar_contains_intron(cigar) and \
                    not self._check_min_contiguous_match(cigar):
                return CIGAR_FAIL
            if self._get_cigar_contains_insertion(cigar) or \
                    self._get_cigar_contains_deletion(cigar):
                return CIGAR_FAIL

        return graded_response

    def _get_cigar_contains_intron(self, cigar):
        for operation, length in cigar:
            if operation == CIGAR_OP_REF_SKIP:
                return True

        return False

    def _get_cigar_contains_insertion(self, cigar):
        for operation, length in cigar:
            if operation == CIGAR_OP_REF_INSERTION:
                return True

        return False

    def _get_cigar_contains_deletion(self, cigar):
        for operation, length in cigar:
            if operation == CIGAR_OP_REF_DELETION:
                return True

        return False


    def _get_cigar_total_match_length(self, cigar):
        total_match = 0

        for operation, length in cigar:
            if operation == CIGAR_OP_MATCH:
                total_match += length

        return total_match


    def _check_min_contiguous_match(self, cigar):
        for operation, length in cigar:
            if operation == CIGAR_OP_MATCH and \
                    length < self.overhang_threshold:
                return False

        return True
