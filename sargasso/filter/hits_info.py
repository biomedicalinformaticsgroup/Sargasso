from sargasso.utlis.factory import Manager


class HitsInfo:
    def __init__(self, data_type, hits):
        self.data_type = data_type
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

    def _is_primary_hit(self, hit):
        return not hit.is_secondary

    def _is_paired_hit(self, hit):
        return hit.is_paired

    def _get_read_length(self, hit):
        return hit.query_length

    def _get_total_length(self, primary_hits):
        total_length = 0
        for hit in primary_hits:
            total_length += self._get_read_length(hit)
        return total_length

    def _get_cigar(self, hit):
        return hit.cigartuples

    def _get_multimaps(self, hits):
        raise NotImplementedError()

    def _get_mismatches(self, hit):
        raise NotImplementedError()

    def _get_alignment_scores(self, hit):
        raise NotImplementedError()


class RnaseqHitsInfo(HitsInfo):

    def _get_multimaps(self, hits):
        return hits[0].get_tag("NH")

    def _get_mismatches(self, hit):
        return hit.get_tag("nM")

    def _get_alignment_scores(self, hit):
        return hit.get_tag("AS")


class ChipseqHitsInfo(HitsInfo):

    def _get_multimaps(self, hits):
        if self._is_paired_hit(hits[0]):
            return len(hits) / 2
        return len(hits)

    def _get_mismatches(self, hit):
        return hit.get_tag("XM")

    def _get_alignment_scores(self, hit):
        return hit.get_tag("AS")


class HitsInfoManager(Manager):
    HITSINFO = {'chipseq': ChipseqHitsInfo,
                'rnaseq': RnaseqHitsInfo}

    @staticmethod
    def get(data_type, hits):
        return HitsInfoManager.HITSINFO[data_type](data_type, hits)