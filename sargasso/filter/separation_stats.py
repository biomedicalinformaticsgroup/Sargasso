class SeparationStats:
    def __init__(self, species_id):
        self.name = "Species {n}".format(n=species_id)
        self.hits_written = 0
        self.reads_written = 0
        self.hits_rejected = 0
        self.reads_rejected = 0
        self.hits_ambiguous = 0
        self.reads_ambiguous = 0

    def accepted_hits(self, hits):
        self.hits_written += len(hits)
        self.reads_written += 1

    def rejected_hits(self, hits):
        self.hits_rejected += len(hits)
        self.reads_rejected += 1

    def ambiguous_hits(self, hits):
        self.hits_ambiguous += len(hits)
        self.reads_ambiguous += 1

    def __str__(self):
        return ("{n}: wrote {f} filtered hits for {fr} reads; {r} hits for " +
                "{rr} reads were rejected outright, and {a} hits for " +
                "{ar} reads were rejected as ambiguous.").format(
            n=self.name,
            f=self.hits_written, fr=self.reads_written,
            r=self.hits_rejected, rr=self.reads_rejected,
            a=self.hits_ambiguous, ar=self.reads_ambiguous)
