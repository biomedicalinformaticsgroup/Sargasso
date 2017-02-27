import sys

from . import samutils


class HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.read_length = None
        self.total_length = None
        self.multimaps = None
        self.max_mismatches = None
        self.min_mismatches = None
        self.cigar = None

    def get_read_length(self):
        if self.read_length is None:
            self.read_length = samutils.get_read_length(self.hits[0])
        return self.read_length

    def get_total_length(self):
        if self.total_length is None:
            self.total_length = samutils.get_total_length(self.hits[0])
        return self.total_length

    def get_multimaps(self):
        if self.multimaps is None:
            self.multimaps = samutils.get_multimaps(self.hits[0])
        return self.multimaps

    def get_max_mismatches(self):
        if self.max_mismatches is None:
            self.max_mismatches = 0
            for hit in self.hits:
                mismatches = samutils.get_mismatches(hit)
                if mismatches > self.max_mismatches:
                    self.max_mismatches = mismatches
        return self.max_mismatches

    def get_min_mismatches(self):
        if self.min_mismatches is None:
            self.min_mismatches = sys.maxsize
            for hit in self.hits:
                mismatches = samutils.get_mismatches(hit)
                if mismatches < self.min_mismatches:
                    self.min_mismatches = mismatches
        return self.min_mismatches

    def get_cigars(self):
        if self.cigar is None:
            self.cigar = {}
            for hit in self.hits:
                self.cigar[hit] = samutils.get_cigar(hit)
        return self.cigar.values()
