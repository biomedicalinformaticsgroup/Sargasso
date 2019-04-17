class HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.primary_hits = self._get_primary_hits()
        self.total_length = self._get_total_length(self.primary_hits)
        self.primary_cigars = [self._get_cigar(h) for h in self.primary_hits]
        self.multimaps = self._get_multimaps(self.hits)
        self.primary_mismatches = self._get_mismatches(self.primary_hits[0])

    def get_total_length(self):
        return self.total_length

    def get_multimaps(self):
        return self.multimaps

    def get_primary_mismatches(self):
        return self.primary_mismatches

    def get_primary_cigars(self):
        return self.primary_cigars

    def _get_primary_hits(self):
        first_hit = None
        for hit in self.hits:
            if self._is_primary_hit(hit):
                if self._is_paired_hit(hit):
                    if first_hit is None:
                        first_hit = hit
                    else:
                        return [first_hit, hit]
                else:
                    return [hit]

        return [first_hit]

    @classmethod
    def _is_primary_hit(cls, hit):
        return not hit.is_secondary

    @classmethod
    def _is_paired_hit(cls, hit):
        return hit.is_paired

    @classmethod
    def _get_read_length(cls, hit):
        return hit.query_length

    @classmethod
    def _get_total_length(cls, primary_hits):
        total_length = 0
        for hit in primary_hits:
            total_length += cls._get_read_length(hit)
        return total_length

    @classmethod
    def _get_cigar(cls, hit):
        return hit.cigartuples

    @classmethod
    def _get_multimaps(cls, hits):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def _get_mismatches(cls, hit):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def _get_alignment_scores(cls, hit):
        raise NotImplementedError('Need to implement in subclass')


class RnaSeqHitsInfo(HitsInfo):
    @classmethod
    def _get_multimaps(cls, hits):
        return hits[0].get_tag("NH")

    @classmethod
    def _get_mismatches(cls, hit):
        return hit.get_tag("nM")

    @classmethod
    def _get_alignment_scores(cls, hit):
        return hit.get_tag("AS")


class DnaSeqHitsInfo(HitsInfo):
    @classmethod
    def _get_multimaps(cls, hits):
        if cls._is_paired_hit(hits[0]):
            return len(hits) / 2
        return len(hits)

    @classmethod
    def _get_mismatches(cls, hit):
        return hit.get_tag("XM")

    @classmethod
    def _get_alignment_scores(cls, hit):
        return hit.get_tag("AS")


class BisulfiteHitsInfo(HitsInfo):
#     Bismark BAM/SAM OUTPUT (default):
#
# (1) QNAME  (seq-ID)
# (2) FLAG   (this flag tries to take the strand a bisulfite read originated from into account (this is different from ordinary DNA alignment flags!))
# (3) RNAME  (chromosome)
# (4) POS    (start position)
# (5) MAPQ   (always 255 for use with Bowtie)
# (6) CIGAR
# (7) RNEXT
# (8) PNEXT
# (9) TLEN
# (10) SEQ
# (11) QUAL   (Phred33 scale)
# (12) NM-tag (edit distance to the reference)
# (13) MD-tag (base-by-base mismatches to the reference (handles indels)
# (14) XM-tag (methylation call string)
# (15) XR-tag (read conversion state for the alignment)
# (16) XG-tag (genome conversion state for the alignment)
# (17) XA/XB-tag (non-bisulfite mismatches) (optional!)
    #todo bisulfite
    pass
