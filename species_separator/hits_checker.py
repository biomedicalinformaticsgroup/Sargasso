ASSIGNED_TO_SPECIES_ONE = 1
ASSIGNED_TO_SPECIES_TWO = 2
ASSIGNED_TO_NEITHER_REJECTED = 3
ASSIGNED_TO_NEITHER_AMBIGUOUS = 4

CIGAR_GOOD = 0
CIGAR_LESS_GOOD = 1
CIGAR_FAIL = 2


class HitsChecker:
    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh, logger):
        self.mismatch_thresh = int(mismatch_thresh)
        self.minmatch_thresh = int(minmatch_thresh)
        self.multimap_thresh = int(multimap_thresh)

        logger.debug(("PARAMS: mismatch - {mism}, minmatch - {minm}, " +
                      "multimap - {mult}").format(
            mism=mismatch_thresh,
            minm=minmatch_thresh,
            mult=multimap_thresh))

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species. For example, in the ultra-conservative
        # strategy, this means:
        # - no multi-mapping
        # - no mismatches
        # - CIGAR string satisfactory

        if hits_info.get_multimaps() > self.multimap_thresh:
            return False

        if hits_info.get_max_mismatches() > self.mismatch_thresh:
            return False

        if self.check_cigar(hits_info) == CIGAR_FAIL:
            return False

        return True

    def assign_hits(self, s1_hits_info, s2_hits_info):
        # TODO: currently making the assumption that min = max mismatches
        # TODO: initially try to get this to behave in the same way as the original code

        violated1 = False
        violated2 = False

        s1_multimaps = s1_hits_info.get_multimaps()
        s2_multimaps = s2_hits_info.get_multimaps()

        if s1_multimaps > self.multimap_thresh:
            violated1 = True
        if s2_multimaps > self.multimap_thresh:
            violated2 = True

        s1_mismatches = s1_hits_info.get_max_mismatches()
        s2_mismatches = s2_hits_info.get_max_mismatches()

        if s1_mismatches > self.mismatch_thresh:
            violated1 = True
        if s2_mismatches > self.mismatch_thresh:
            violated2 = True

        s1_cigar_check = self.check_cigar(s1_hits_info)
        s2_cigar_check = self.check_cigar(s2_hits_info)

        if s1_cigar_check == CIGAR_FAIL:
            violated1 = True
        if s2_cigar_check == CIGAR_FAIL:
            violated2 = True

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
        elif violated2 and violated1:
            return ASSIGNED_TO_NEITHER_REJECTED

        return ASSIGNED_TO_NEITHER_AMBIGUOUS

    def retrieve_full_number(self, cigar):
        for i in range(len(cigar)):
            if len(cigar) - (i + 1) > 0:
                if not cigar[len(cigar) - (i + 1)].isdigit():
                    return cigar[(len(cigar) - i):len(cigar)]
            else:
                return cigar[0:len(cigar)]
        return 0  # should never reach here, but just in case

    # extracts the quantity of a certain type of base; e.g M, S, N etc
    def extract_base_quantity(self, cigar, type):
        if not type in cigar:
            return 0
        else:
            total = 0
            for i in range(len(cigar)):
                if cigar[i] == type:
                    total += int(self.retrieve_full_number(cigar[0:i]))
            return total

    def check_min_match(self, cigar):
        if not "M" in cigar:
            return False
        else:
            quantity = 0
            for i in range(len(cigar)):
                if cigar[i] == "M":
                    quantity = int(self.retrieve_full_number(cigar[0:i]))
                    if quantity < 5:
                        return False
            return True

    def check_cigar(self, hits_info):
        length = hits_info.get_length()
        min_match = length - self.minmatch_thresh

        # when allowing other params, this is no longer a t/f scenario;
        # e.g. when clipping is allowed a cigar without clipping should score
        # better than one with even though both are allowed
        graded_response = CIGAR_GOOD

        for cigar in hits_info.get_cigars():
            # DO NOT REMOVE:
            #for c in ["I", "D", "S", "H", "P", "X"]: # ENABLE FOR CONSERVATIVE
            #    if c in cigar:
            #        return 2

            # ENABLE FOR MINMATCH THRESHOLDING:
            num_matches = self.extract_base_quantity(cigar, "M")
            if num_matches < min_match:
                return CIGAR_FAIL
            elif num_matches < length:
                graded_response = CIGAR_LESS_GOOD

            if "N" in cigar:
                if not self.check_min_match(cigar):
                    return CIGAR_FAIL

        return graded_response
