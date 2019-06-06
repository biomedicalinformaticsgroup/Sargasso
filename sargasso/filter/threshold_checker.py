import json


class ThresholdChecker(object):
    CIGAR_GOOD = 0
    CIGAR_LESS_GOOD = 1
    CIGAR_FAIL = 2

    CIGAR_OP_MATCH = 0  # From pysam
    CIGAR_OP_REF_INSERTION = 1  # From pysam
    CIGAR_OP_REF_DELETION = 2  # From pysam
    CIGAR_OP_REF_SKIP = 3  # From pysam

    def __init__(self, hits_manager, mismatch_thresh, minmatch_thresh, multimap_thresh):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species.
        # return True if valid

        self.species_id = hits_manager.species_id
        hits_info = hits_manager.hits_info

        self._check_multimap_thresh(hits_info, multimap_thresh)
        self._check_mismatch_thresh(hits_info, mismatch_thresh)
        self._check_cigar_thresh(hits_info, minmatch_thresh)
        self._check_violated()

    def _check_multimap_thresh(self, hits_info, multimap_thresh):
        self.violated_multimaps = False
        self.multimaps = hits_info.get_multimaps()
        if self.multimaps > multimap_thresh:
            self.violated_multimaps = True

    def _check_mismatch_thresh(self, hits_info, mismatch_thresh):
        self.violated_mismatches = False
        self.mismatches = hits_info.get_primary_mismatches()
        if self.mismatches > round(mismatch_thresh * hits_info.get_total_length()):
            self.violated_mismatches = True

    def _check_cigar_thresh(self, hits_info, minmatch_thresh):
        self.violated_cigar_check = False
        self.cigar_check = self._check_cigars(hits_info, minmatch_thresh)
        if self.cigar_check == self.CIGAR_FAIL:
            self.violated_cigar_check = True

    def _check_violated(self):
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

    def to_debug_string(self):
        dic = self._to_debug_dict()
        return json.dumps(dic)

    def _to_debug_dict(self):
        return {'species_id': self.species_id, 'violated': self.violated,
                'multimaps': self.multimaps, 'mismatches': self.mismatches, 'cigar_check': self.cigar_check,
                'violated_multimaps': self.violated_multimaps,
                'violated_mismatches': self.violated_mismatches,
                'violated_cigar_check': self.violated_cigar_check
                }


class RnaThresholdChecker(ThresholdChecker):
    pass


class DnaThresholdChecker(ThresholdChecker):
    pass


class BisulfiteThresholdChecker(ThresholdChecker):
    def __init__(self, hits_manager, mismatch_thresh, minmatch_thresh, multimap_thresh):
        super(self.__class__, self).__init__(hits_manager, mismatch_thresh, minmatch_thresh, multimap_thresh)

        self._check_ambig(hits_manager.hits_info)
        self.violated = self.violated or self.violated_ambig

    def _check_ambig(self, hits_info):
        self.is_ambig = hits_info.get_is_ambig_hit()
        self.violated_ambig = self.is_ambig

    def _to_debug_dict(self):
        debug_list = super(self.__class__, self)._to_debug_dict()
        debug_list['violated_ambig'] = self.violated_ambig
        debug_list['is_ambig'] = self.is_ambig
        return debug_list
