ASSIGNED_TO_SPECIES_ONE = 1
ASSIGNED_TO_SPECIES_TWO = 2
ASSIGNED_TO_NEITHER_REJECTED = 3
ASSIGNED_TO_NEITHER_AMBIGUOUS = 4

CIGAR_GOOD = 0
CIGAR_LESS_GOOD = 1
CIGAR_FAIL = 2

REJECT_EDITS_OPERATIONS = ["I", "D", "S", "H", "P", "X"]


class HitsChecker:
    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh,
                 reject_multimaps, reject_edits, logger):
        self.mismatch_thresh = mismatch_thresh
        self.minmatch_thresh = minmatch_thresh
        self.multimap_thresh = multimap_thresh
        self._assign_hits = self._assign_hits_reject_multimaps \
            if reject_multimaps else self._assign_hits_standard
        self._check_cigar = self._check_cigar_reject_edits \
            if reject_edits else self._check_cigar_standard

        logger.debug(("PARAMS: mismatch - {mism}, minmatch - {minm}, " +
                      "multimap - {mult}").format(
            mism=mismatch_thresh, minm=minmatch_thresh, mult=multimap_thresh))

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        return hits_info.get_multimaps() <= self.multimap_thresh and \
            hits_info.get_max_mismatches() <= self.mismatch_thresh and \
            self._check_cigars(hits_info) != CIGAR_FAIL

    def assign_hits(self, s1_hits_info, s2_hits_info):
        # TODO: currently making the assumption that min = max mismatches
        violated1, s1_multimaps, s1_mismatches, s1_cigar_check = \
            self._check_thresholds(s1_hits_info)
        violated2, s2_multimaps, s2_mismatches, s2_cigar_check = \
            self._check_thresholds(s2_hits_info)

        return self._assign_hits(
            violated1, s1_multimaps, s1_mismatches, s1_cigar_check,
            violated2, s2_multimaps, s2_mismatches, s2_cigar_check)

    def _assign_hits_standard(
            self, violated1, s1_multimaps, s1_mismatches, s1_cigar_check,
            violated2, s2_multimaps, s2_mismatches, s2_cigar_check):

        if violated1 and violated2:
            return ASSIGNED_TO_NEITHER_REJECTED

        if s1_mismatches < s2_mismatches and not violated1:
            return ASSIGNED_TO_SPECIES_ONE
        elif s2_mismatches < s1_mismatches and not violated2:
            return ASSIGNED_TO_SPECIES_TWO
        else:
            if s1_cigar_check < s2_cigar_check and not violated1:
                return ASSIGNED_TO_SPECIES_ONE
            elif s2_cigar_check < s1_cigar_check and not violated2:
                return ASSIGNED_TO_SPECIES_TWO
            else:
                if s1_multimaps < s2_multimaps and not violated1:
                    return ASSIGNED_TO_SPECIES_ONE
                elif s2_multimaps < s1_multimaps and not violated2:
                    return ASSIGNED_TO_SPECIES_TWO

        if not violated1 and violated2:
            return ASSIGNED_TO_SPECIES_ONE
        elif not violated2 and violated1:
            return ASSIGNED_TO_SPECIES_TWO

        return ASSIGNED_TO_NEITHER_AMBIGUOUS

    def _assign_hits_reject_multimaps(
            self, violated1, s1_multimaps, s1_mismatches, s1_cigar_check,
            violated2, s2_multimaps, s2_mismatches, s2_cigar_check):

        if s1_multimaps > 1 or s2_multimaps > 1:
            return ASSIGNED_TO_NEITHER_REJECTED

        return self._assign_hits_standard(
            violated1, s1_multimaps, s1_mismatches, s1_cigar_check,
            violated2, s2_multimaps, s2_mismatches, s2_cigar_check)

    def _check_thresholds(self, hits_info):
        violated = False

        multimaps = hits_info.get_multimaps()
        if multimaps > self.multimap_thresh:
            violated = True

        mismatches = hits_info.get_max_mismatches()
        if mismatches > self.mismatch_thresh:
            violated = True

        cigar_check = self._check_cigars(hits_info)
        if cigar_check == CIGAR_FAIL:
            violated = True

        return (violated, multimaps, mismatches, cigar_check)

    def _check_cigars(self, hits_info):
        length = hits_info.get_length()
        min_match = length - self.minmatch_thresh

        # when allowing other params, this is no longer a t/f scenario;
        # e.g. when clipping is allowed a cigar without clipping should score
        # better than one with even though both are allowed
        graded_response = CIGAR_GOOD

        for cigar in hits_info.get_cigars():
            graded_response = self._check_cigar(
                cigar, graded_response, min_match, length)
            if graded_response == CIGAR_FAIL:
                return graded_response

        return graded_response

    def _check_cigar_standard(
            self, cigar, current_response, min_match, length):

        graded_response = current_response

        num_matches = self._get_total_cigar_op_length(cigar, "M")
        if num_matches < min_match:
            return CIGAR_FAIL
        elif num_matches < length:
            graded_response = CIGAR_LESS_GOOD

        if "N" in cigar:
            if not self._check_min_contiguous_match(cigar):
                return CIGAR_FAIL

        return graded_response

    def _check_cigar_reject_edits(
            self, cigar, current_response, min_match, length):

        for operation in REJECT_EDITS_OPERATIONS:
            if operation in cigar:
                return CIGAR_FAIL

        return self._check_cigar_standard(
            cigar, current_response, min_match, length)

    # extracts the quantity of a certain type of base; e.g M, S, N etc
    def _get_total_cigar_op_length(self, cigar, op_type):
        total = 0
        if op_type in cigar:
            for i, char in enumerate(cigar):
                if char == op_type:
                    total += self._get_length_of_cigar_op(cigar, i)

        return total

    def _check_min_contiguous_match(self, cigar):
        if not "M" in cigar:
            return False

        for i, char in enumerate(cigar):
            if char == "M":
                if self._get_length_of_cigar_op(cigar, i) < 5:
                    return False
        return True

    def _get_length_of_cigar_op(self, cigar, index):
        for i in range(index - 1):
            if not cigar[index - i - 1].isdigit():
                return int(cigar[(index - i):index])
        return int(cigar[:index])