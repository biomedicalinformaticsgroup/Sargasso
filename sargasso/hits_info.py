import sys

from . import samutils


class HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.primary_hits = self._get_primary_hits()
        self.total_length = samutils.get_total_length(self.primary_hits)
        self.multimaps = samutils.get_multimaps(self.hits[0])
        self.primary_mismatches = samutils.get_mismatches(self.primary_hits[0])
        self.primary_cigars = [samutils.get_cigar(h) for h in self.primary_hits]

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
            if samutils.is_primary(hit):
                if hit.is_paired:
                    if first_hit is None:
                        first_hit = hit
                    else:
                        return [first_hit, hit]
                else:
                    return [hit]

        return [first_hit]
