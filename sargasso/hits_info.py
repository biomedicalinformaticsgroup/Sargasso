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
            self.primary_cigars = {}
            for hit in self.get_primary_hits():
                self.primary_cigars[hit] = samutils.get_cigar(hit)
        return self.primary_cigars.values()

    def get_primary_hits(self):
        if self.primary_hits is None:
            self.primary_hits = []
            is_paried = self.hits[0].is_paired
            for hit in self.hits:
                if samutils.is_primary(hit):
                    self.primary_hits.append(hit)
                    #to save some interation if we find all the primary hits
                    #before reaching the end of the list
                    if not is_paried or len(self.primary_hits) == 2:
                        return self.primary_hits
        return self.primary_hits



