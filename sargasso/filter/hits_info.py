class HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.primary_hits = self._get_primary_hits()
        self.total_length = self._get_total_length(self.primary_hits)
        self.primary_cigars = [self._get_cigar(h) for h in self.primary_hits]
        self.multimaps = self._get_num_multimaps(self.hits)
        self.primary_mismatches = self._get_mismatches(self.primary_hits)

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
    def _get_num_multimaps(cls, hits):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def _get_mismatches(cls, hits):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def _get_alignment_scores(cls, hit):
        return hit.get_tag("AS")


class RnaSeqHitsInfo(HitsInfo):
    @classmethod
    def _get_num_multimaps(cls, hits):
        return hits[0].get_tag("NH")

    @classmethod
    def _get_mismatches(cls, hits):
        return hits[0].get_tag("nM")




class DnaSeqHitsInfo(HitsInfo):
    @classmethod
    def _get_num_multimaps(cls, hits):
        if cls._is_paired_hit(hits[0]):
            return len(hits) / 2
        return len(hits)

    @classmethod
    def _get_mismatches(cls, hits):
        return hits[0].get_tag("XM")



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

    def __init__(self, hits):
        HitsInfo.__init__(self,hits)
        # @xintodo The reads in the pair can have a different NM field.
        # self.primary_edit_distance = self._get_edit_distance(self.primary_hits)
        self.is_ambig_hit=self._is_ambig_hit(self.primary_hits[0])


    #todo bisulfite
    # This should always be 1 as bowtie2 return only the best hit as configured by bismark
    @classmethod
    def _get_num_multimaps(cls, hits):
        if cls._is_paired_hit(hits[0]):
            return len(hits) / 2
        return len(hits)

    ## bismark only return one (pair of) hit(s) for each read
    ## we want to know if this hit is A alignment from multi-maps (in bismark ambig file, contains AS: field)
    ## or a singleton alignment (in bismark output bam file, does NOT contain AS:field)
    @classmethod
    def _is_ambig_hit(cls, hit):
        return hit.has_tag("AS")


    #NM:i:<N> The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and deletions) needed to transform the read string into the reference string. Only present if SAM record is for an aligned read.
    # @classmethod
    # def _get_edit_distance(cls, hits):
    #     if cls._is_paired_hit(hits[0]):
    #         sum([hit.get_tag("NM") for hit in hits])
    #     return hits[0].get_tag("NM")

    # --non_bs_mm Optionally outputs an extra column specifying the number of non-bisulfite
    #               mismatches a read during the alignment step. This option is only available for SAM
    #               format. In Bowtie 2 context, this value is just the number of actual non-bisulfite
    #               mismatches and ignores potential insertions or deletions. The format for single-end
    #               reads and read 1 of paired-end reads is 'XA:Z:number of mismatches' and
    #               'XB:Z:number of mismatches' for read 2 of paired-end reads.
    @classmethod
    def _get_mismatches(cls, hits):
        return max([cls._get_mismatches_by_hit(hit) for hit in hits])


    @classmethod
    def _get_mismatches_by_hit(cls,hit):
        if cls._is_ambig_hit(hit):
            return cls._get_mismatches_ambig(hit)

        if hit.is_read1 or not cls._is_paired_hit(hit):
            return cls._get_mismatches_read1(hit)

        return cls._get_mismatches_read2(hit)

    @classmethod
    def _get_mismatches_read1(cls,hit):
        return int(hit.get_tag("XA"))

    @classmethod
    def _get_mismatches_read2(cls,hit):
        return int(hit.get_tag("XB"))

    @classmethod
    def _get_mismatches_ambig(cls,hit):
        return int(hit.get_tag("XM"))

    def get_is_ambig_hit(self):
        return self.is_ambig_hit
