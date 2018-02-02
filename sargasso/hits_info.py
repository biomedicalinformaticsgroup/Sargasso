import sys

from . import samutils


class HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.total_length = None
        self.multimaps = None
        self.primary_mismatches = None
        self.primary_cigars = None
        self.primary_hits = None

    def get_total_length(self):
        if self.total_length is None:
            self.total_length = samutils.get_total_length(self.get_primary_hits())
        return self.total_length

    def get_multimaps(self):
        if self.multimaps is None:
            self.multimaps = samutils.get_multimaps(self.hits[0])
        return self.multimaps

    def get_primary_mismatches(self):
        if self.primary_mismatches is None:
            primary_hit = self.get_primary_hits()[0]
            self.primary_mismatches=samutils.get_mismatches(primary_hit)
        return self.primary_mismatches


    def get_primary_cigars(self):
            if self.primary_cigars is None:
                self.primary_cigars = [samutils.get_cigar(h) for h in self.get_primary_hits()]
            return self.primary_cigars

    def get_primary_hits(self):
        if self.primary_hits is None:
            first_hit = None
            for hit in self.hits:
                if samutils.is_primary(hit):
                    if hit.is_paired:
                        if first_hit is None:
                            first_hit = hit
                        else:
                            self.primary_hits = [first_hit, hit]
                            return self.primary_hits
                    else:
                        self.primary_hits = [hit]
                        return self.primary_hits
            self.primary_hits = [first_hit]
        return self.primary_hits



