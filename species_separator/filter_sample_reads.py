#!/usr/bin/env python

"""Usage:
    filter_mapped_hits [--log-level=<log-level>] <species-one> <species-one-input-bam> <species-one-output-bam> <species-two> <species-two-input-bam> <species-two-output-bam> <mismatch-threshold> <minmatch-threshold> <multimap-threshold> <alignment-score-threshold>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>               Name of first species.
<species-one-input-bam>     BAM file containing read hits against first species' genome.
<species-one-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<species-two>               Name of first species.
<species-two-input-bam>     BAM file containing read hits against first species' genome.
<species-two-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<mismatch-threshold>	    Maximum number of mismatches allowed during filtering
<minmatch-threshold>        Minimum number of read bases that must be perfectly matched
<multimap-threshold>        Maximum number of multiple mappings allowed during filtering
<alignment-score-threshold> Minimum alignment score allowed during filtering

TODO: what does this script do...
"""

import docopt
import schema
import sys

from . import options as opt
from . import samutils
from .__init__ import __version__

SPECIES_ONE = "<species-one>"
SPECIES_ONE_INPUT_BAM = "<species-one-input-bam>"
SPECIES_ONE_OUTPUT_BAM = "<species-one-output-bam>"
SPECIES_TWO = "<species-two>"
SPECIES_TWO_INPUT_BAM = "<species-two-input-bam>"
SPECIES_TWO_OUTPUT_BAM = "<species-two-output-bam>"
MISMATCH_THRESHOLD = "<mismatch-threshold>"
MINMATCH_THRESHOLD = "<minmatch-threshold>"
MULTIMAP_THRESHOLD = "<multimap-threshold>"
ALIGNMENT_SCORE_THRESHOLD = "<alignment-score-threshold>"

ASSIGNED_TO_SPECIES_ONE = 1
ASSIGNED_TO_SPECIES_TWO = 2
ASSIGNED_TO_NEITHER_REJECTED = 3
ASSIGNED_TO_NEITHER_AMBIGUOUS = 4


class _SeparationStats:
    def __init__(self, name):
        self.name = name
        self.hits_written = 0
        self.reads_written = 0
        self.hits_rejected = 0
        self.reads_rejected = 0
        self.hits_ambiguous = 0
        self.reads_ambiguous = 0

    def filtered_hits(self, hits):
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


def _validate_command_line_options(options):
    try:
        opt.validate_log_level(options)

        opt.validate_file_option(
            options[SPECIES_ONE_INPUT_BAM],
            "Could not find input BAM file for species 1")
        opt.validate_file_option(
            options[SPECIES_TWO_INPUT_BAM],
            "Could not find input BAM file for species 2")
    except schema.SchemaError as exc:
        exit(exc.code)


class _HitsInfo:
    def __init__(self, hits):
        self.hits = hits
        self.multimaps = None
        self.max_mismatches = None
        self.min_mismatches = None
        self.alignment_scores = None
        self.cigar = None

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

    def get_min_alignment_score(self):
        if self.alignment_scores is None:
            self.alignment_scores = sys.maxsize
            for hit in self.hits:
                score = samutils.get_alignment_scores(hit)
                if score < self.alignment_scores:
                    self.alignment_scores = score
        return self.alignment_scores

    def get_cigars(self):
        if self.cigar is None:
            self.cigar = {}
            for hit in self.hits:
                self.cigar[hit] = samutils.get_cigar(hit)
        return self.cigar.values()


class _HitsChecker:
    def __init__(self, mismatch_thresh, minmatch_thresh, multimap_thresh, alignment_score_thresh):
        self.mismatch_thresh = mismatch_thresh
        self.minmatch_thresh = minmatch_thresh
        self.multimap_thresh = multimap_thresh
        self.alignment_score_thresh = alignment_score_thresh
        #print "PARAMS: mismatch - "+str(mismatch_thresh)+", minmatch - "+str(minmatch_thresh)+", multimap - "+str(multimap_thresh)+", alignment score - "+str(alignment_score_thresh)

    def check_hits(self, hits_info):
        # check that the hits for a read are - in themselves - satisfactory to
        # be assigned to a species. For example, in the ultra-conservative
        # strategy, this means:
        # - no multi-mapping
        # - no mismatches
        # - CIGAR string satisfactory

        #if hits_info.get_multimaps() > 1:
        if hits_info.get_multimaps() > self.multimap_thresh: #MULTIMAP_THRESHOLD:
            return False

        #if hits_info.get_max_mismatches() > 0:
        if hits_info.get_max_mismatches() > self.mismatch_thresh: #MISMATCH_THRESHOLD:
            return False

        # ENABLE FOR AS FILTERING
        if self.check_alignment_score(hits_info) == 2:
            return False

        #if not self.check_cigar(hits_info):
        if self.check_cigar(hits_info) == 2:
            return False

        return True

    def assign_hits(self, s1_hits_info, s2_hits_info):
        # TODO: currently making the assumption that min = max mismatches
        # TODO: initially try to get this to behave in the same way as the original code

        s1_multimaps = s1_hits_info.get_multimaps()
        s2_multimaps = s2_hits_info.get_multimaps()

        mmThresh = self.multimap_thresh #MULTIMAP_THRESHOLD

        if s1_multimaps > mmThresh or s2_multimaps > mmThresh:
            return ASSIGNED_TO_NEITHER_REJECTED

        #print "NEW Param Check: Mismatches - "+str(options[MISMATCH_THRESHOLD])+", Minmatch - "+str(options[MINMATCH_THRESHOLD])+", Multimaps - "+str(options[MULTIMAP_THRESHOLD])+", AS - "+str(options[ALIGNMENT_SCORE_THRESHOLD])

	s1_AS_check = self.check_alignment_score(s1_hits_info)
        s2_AS_check = self.check_alignment_score(s2_hits_info)
 
        if s1_AS_check==0 and s2_AS_check==2:
            return ASSIGNED_TO_SPECIES_ONE
        elif s2_AS_check==0 and s1_AS_check==2:
            return ASSIGNED_TO_SPECIES_TWO

        s1_mismatches = s1_hits_info.get_max_mismatches()
        s2_mismatches = s2_hits_info.get_max_mismatches()

        if s1_mismatches < s2_mismatches:
            return ASSIGNED_TO_SPECIES_ONE
        elif s2_mismatches < s1_mismatches:
            return ASSIGNED_TO_SPECIES_TWO
        else:
            s1_cigar_check = self.check_cigar(s1_hits_info)
            s2_cigar_check = self.check_cigar(s2_hits_info)

            #if s1_cigar_check and not s2_cigar_check:
            if s1_cigar_check<s2_cigar_check:
                return ASSIGNED_TO_SPECIES_ONE
            #elif not s1_cigar_check and s2_cigar_check:
            elif s2_cigar_check<s1_cigar_check:
                return ASSIGNED_TO_SPECIES_TWO
            
            #s1_AS_check = self.check_alignment_score(s1_hits_info)
            #s2_AS_check = self.check_alignment_score(s2_hits_info)
 
            #if s1_AS_check==0 and s2_AS_check==2:
            #    return ASSIGNED_TO_SPECIES_ONE
            #elif s2_AS_check==0 and s1_AS_check==2:
            #    return ASSIGNED_TO_SPECIES_TWO

        #print "ONE: "+str(s1_cigar_check)+", TWO: "+str(s2_cigar_check)
        return ASSIGNED_TO_NEITHER_AMBIGUOUS

    def retreive_full_number(self,cigar):
        for i in range(len(cigar)):
            if len(cigar)-(i+1) > 0:
                #print "CIGAR: "+cigar+", i: "+str(i)+", length: "+str(len(cigar))
                if not cigar[len(cigar)-(i+1)].isdigit():
                    #print "CIGAR: "+cigar+", i: "+str(i)+", FINISHED: "+cigar[(len(cigar)-i):len(cigar)]
                    return cigar[(len(cigar)-i):len(cigar)]
            else:
                #print "CIGAR: "+cigar+", FINISHED: "+cigar[0:len(cigar)]
                return cigar[0:len(cigar)]
        return 0 # should never reach here, but just in case

    # extracts the quantity of a certain type of base; e.g M, S, N etc
    def extract_base_quantity(self,cigar,type):
        if not type in cigar:
            return 0
        else:
            total = 0
            for i in range(len(cigar)):
                if cigar[i] == type:
                    total += int(self.retreive_full_number(cigar[0:i]))
                    #print "CIGAR: "+cigar+", Passing: "+cigar[0:i]+", Total: "+str(total)
            return total

    def check_min_match(self,cigar):
        if not "M" in cigar:
            return False
        else:
            quantity = 0
            for i in range(len(cigar)):
                if cigar[i] == "M":
                    quantity = int(self.retreive_full_number(cigar[0:i]))
                    if quantity < 5:
                        return False
                    #print "CIGAR: "+cigar+", Passing: "+cigar[0:i]+", Total: "+str(total)
            return True

    def check_cigar(self, hits_info):
        gradedResponse = 0# when allowing other params, this is no longer a t/f scenario; e.g when clipping is allowed a cigar without clipping should score better than one with even though both are allowed
        # Grades: 2 = fail, 1 = less good, 0 = good
        for cigar in hits_info.get_cigars():
            for c in ["I", "D", "S", "H", "P", "X"]: # ENABLE FOR CONSERVATIVE
            #for c in ["I", "D", "H", "P", "X"]: # ENABLE FOR SOFT CLIPPING
                if c in cigar:
                    return 2

            # ENABLE FOR SOFT CLIPPING
            #allowedClipping = 5 # FO: allowing soft clipping
            #noS = self.extract_base_quantity(cigar,"S")
            #if noS > allowedClipping:
            #    return 2
            #elif noS > 0:
            #    gradedResponse = 1
            
            # ENABLE FOR MINMATCH THRESHOLDING:
            #minMatch = self.minmatch_thresh #MINMATCH_THRESHOLD #40
            #noM = self.extract_base_quantity(cigar,"M")
            #if noM < minMatch:
            #    return 2
            #elif noM < 50:
            #    gradedResponse = 1

            if "N" in cigar:
                #if "S" in cigar: # FILTER OPTIMISING; this allows some soft clipping
                if not self.check_min_match(cigar):
                    return 2
                
                #    segments = cigar.split("S")
                #    if not "M" in segments[0]: # if S only at the start
                #        if int(segments[0]) > allowedClipping:
                #            return 2
                #    elif len(segments)<3:
		#	#print "SEGMENTS: "+str(segments[0])+" "+str(segments[1])
                #        #if int(segments[1]) > allowedClipping:
                #        print "SEGS: "+str(segments)
                #        if int(segments[0].split("M")[2]) > allowedClipping:
                #            return 2
                #    if len(segments)>2: 
                #        if int(segments[1].split("M")[2]) > allowedClipping:
                #            return 2
                #    gradedResponse = 1
                
                #else:
                #segments = cigar.split("M")
                #noN = self.extract_base_quantity(cigar,"N")
		#print "CIGAR: "+cigar+", SEGMENTS: "+str(segments)+", NoN: "+str(noN)
                #if int(segments[0]) < 5 or int(segments[1].split("N")[1]) < 5:
                #if noN <
                #    return 2
        return gradedResponse

    # check whether a read is mapped over a splice junction and if not, whether its score is above threshold
    def check_alignment_score(self,hits_info):
        hit_index = -1
        score_threshold = ALIGNMENT_SCORE_THRESHOLD #100
        cigars = hits_info.get_cigars()
        #for score in hits_info.get_min_alignment_score():
        #print "CIGARS: "+str(cigars)
	#quit()
        score = hits_info.get_min_alignment_score()
        for cig in cigars:
            hit_index += 1
            print "AS: "+str(score)+", CIGAR: "+cigars[hit_index]
            if "N" in cigars[hit_index]:
                return 1 # Not concerning ourselves with these reads for this filtering as we know this filtering doesn't apply to them; as such they are subject to regular filtering
        if score < score_threshold:
            return 2
        else:
            return 0

def _hits_generator(all_hits):
    last_hit_name = None
    current_hits = []

    for hit in samutils.all_hits(all_hits):
        if last_hit_name is None:
            last_hit_name = hit.query_name

        if hit.query_name < last_hit_name:
            # TODO: throw an exception if hits are out of read name order, and
            # handle this gracefully
            pass
        elif hit.query_name == last_hit_name:
            current_hits.append(hit)
        else:
            yield current_hits
            last_hit_name = hit.query_name
            current_hits = [hit]

    yield current_hits


def _write_hits(hits_for_read, output_bam):
    for hit in hits_for_read:
        output_bam.write(hit)
        pass


def _write_remaining_hits(hits_checker, hits_for_reads, stats, output_bam):
    try:
        while True:
            hits_for_read = hits_for_reads.next()
            if hits_checker.check_hits(_HitsInfo(hits_for_read)):
                stats.filtered_hits(hits_for_read)
                _write_hits(hits_for_read, output_bam)
            else:
                stats.rejected_hits(hits_for_read)
    except StopIteration:
        pass


def _compare_and_write_hits(s1_hits_for_read, s2_hits_for_read,
                            s1_output_bam, s2_output_bam,
                            s1_stats, s2_stats, hits_checker):
    # Compare the hits for a particular read in each species and decide whether
    # the read can be assigned to one species or another, or if it must be
    # rejected as ambiguous
    s1_hits_info = _HitsInfo(s1_hits_for_read)
    s2_hits_info = _HitsInfo(s2_hits_for_read)

    assignment = hits_checker.assign_hits(s1_hits_info, s2_hits_info)

    if assignment == ASSIGNED_TO_SPECIES_ONE:
        s2_stats.rejected_hits(s2_hits_for_read)
        if hits_checker.check_hits(s1_hits_info):
            s1_stats.filtered_hits(s1_hits_for_read)
            _write_hits(s1_hits_for_read, s1_output_bam)
        else:
            s1_stats.rejected_hits(s1_hits_for_read)
    elif assignment == ASSIGNED_TO_SPECIES_TWO:
        s1_stats.rejected_hits(s1_hits_for_read)
        if hits_checker.check_hits(s2_hits_info):
            s2_stats.filtered_hits(s2_hits_for_read)
            _write_hits(s2_hits_for_read, s2_output_bam)
        else:
            s2_stats.rejected_hits(s2_hits_for_read)
    elif assignment == ASSIGNED_TO_NEITHER_REJECTED:
        s1_stats.rejected_hits(s1_hits_for_read)
        s2_stats.rejected_hits(s2_hits_for_read)
    else:
        s1_stats.ambiguous_hits(s1_hits_for_read)
        s2_stats.ambiguous_hits(s2_hits_for_read)


def _filter_sample_reads(logger, options):
    logger.info("Starting species separation.")

    s1_stats = _SeparationStats("Species 1")
    s2_stats = _SeparationStats("Species 2")

    s1_hits = samutils.open_samfile_for_read(options[SPECIES_ONE_INPUT_BAM])
    s2_hits = samutils.open_samfile_for_read(options[SPECIES_TWO_INPUT_BAM])

    s1_hits_for_reads = _hits_generator(s1_hits)
    s2_hits_for_reads = _hits_generator(s2_hits)

    s1_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_ONE_OUTPUT_BAM], s1_hits)
    s2_output_bam = samutils.open_samfile_for_write(
        options[SPECIES_TWO_OUTPUT_BAM], s2_hits)

    s1_hits_for_read = None
    s2_hits_for_read = None

    s1_count = 0
    s2_count = 0

    hits_checker = _HitsChecker(options[MISMATCH_THRESHOLD],options[MINMATCH_THRESHOLD],options[MULTIMAP_THRESHOLD],options[ALIGNMENT_SCORE_THRESHOLD])

    while True:
    # 1. Attempt to read all hits for the first/next read in each species'
    # input BAM file. Go to (2).
        if s1_hits_for_read is None:
            try:
                s1_hits_for_read = s1_hits_for_reads.next()
                s1_count += 1
                if s1_count % 1000000 == 0:
                    logger.debug("Read " + str(s1_count) + " reads from species 1")
            except StopIteration:
    # 2. If no more reads can be read for species 1:
    #       - all remaining reads in the input file for species 2 can be written
    #       to the output file for species 2, or discarded as inherently
    #       potentially ambiguous.
    #       - END.
    #    Else go to (3).
                _write_remaining_hits(hits_checker, s2_hits_for_reads, s2_stats, s2_output_bam)
                break

        if s2_hits_for_read is None:
            try:
                s2_hits_for_read = s2_hits_for_reads.next()
                s2_count += 1
                if s2_count % 1000000 == 0:
                    logger.debug("Read " + str(s2_count) + " reads from species 2")
            except StopIteration:
    # 3. If no more reads can be read for species 2:
    #       - all remaining reads in the input file for species 1 can be
    #       written to the output file for species 1, or discarded as
    #       inherently potentially ambiguous.
    #       - END.
    #   Else go to (4).
                _write_remaining_hits(hits_checker, s1_hits_for_reads, s1_stats, s1_output_bam)
                break

    # 4. We have the hits for a read for each species. Examine the names of the
    # reads in species 1 and species 2. Go to (5).
        s1_read_name = s1_hits_for_read[0].query_name
        s2_read_name = s2_hits_for_read[0].query_name

    # 5. If the read names are the same:
        if s1_read_name == s2_read_name:
    #       - compare the hits for each species to determine which species to
    #       assign the read to.
    #       - write the hits to the output file for the correct species (or
    #       neither, if the read cannot be unambiguously assigned).
    #       - Go to (1).
    #   Else go to (6).
            _compare_and_write_hits(s1_hits_for_read, s2_hits_for_read,
                                    s1_output_bam, s2_output_bam,
                                    s1_stats, s2_stats, hits_checker)
            s1_hits_for_read = None
            s2_hits_for_read = None
    # 6. If species 1 read name < species 2 read name, there are no hits for
    # the species 1 read in species 2:
        elif s1_read_name < s2_read_name:
    #       - write the hits for this read to the output file for species 1, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 1. If no more
    #       reads can be read for species 1, go to (2). Else go to (4).
    #   Else go to (7)
            if hits_checker.check_hits(_HitsInfo(s1_hits_for_read)):
                s1_stats.filtered_hits(s1_hits_for_read)
                _write_hits(s1_hits_for_read, s1_output_bam)
            else:
                s1_stats.rejected_hits(s1_hits_for_read)
            s1_hits_for_read = None
    # 7. If species 2 read name < species 1 read name, there are no hits for
    # the species 2 read in species 1:
        else:
    #       - write the hits for this read to the output file for species 2, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 2. If no more
    #       reads can be read for species 2, go to (3). Else go to (4).
            if hits_checker.check_hits(_HitsInfo(s2_hits_for_read)):
                s2_stats.filtered_hits(s2_hits_for_read)
                _write_hits(s2_hits_for_read, s2_output_bam)
            else:
                s2_stats.rejected_hits(s2_hits_for_read)
            s2_hits_for_read = None
    # n.b. throughout, we should always keep a record of the last read name
    # read for a species; when hits for the next read are read for that
    # species, we should check the reads are in name order, and exit with a
    # failure code if not.

    logger.info(s1_stats)
    logger.info(s2_stats)

    write_stats(s1_stats,s2_stats,options[SPECIES_ONE_OUTPUT_BAM])

# write filter stats to table in file
def write_stats(s1,s2,outBam):
    #out = "Filtered-Hits-S1\tFiltered-Reads-S1\tRejected-Hits-S1\tRejected-Reads-S1\tAmbiguous-Hits-S1\tAmbiguous-Reads-S1\tFiltered-Hits-S2\tFiltered-Reads-S2\tRejected-Hits-S2\tRejected-Reads-S2\tAmbiguous-Hits-S2\tAmbiguous-Reads-S2\n"
    out = str(s1.hits_written)+"\t"+str(s1.reads_written)+"\t"+str(s1.hits_rejected)+"\t"+str(s1.reads_rejected)+"\t"+str(s1.hits_ambiguous)+"\t"+str(s1.reads_ambiguous)+"\t"+str(s2.hits_written)+"\t"+str(s2.reads_written)+"\t"+str(s2.hits_rejected)+"\t"+str(s2.reads_rejected)+"\t"+str(s2.hits_ambiguous)+"\t"+str(s2.reads_ambiguous)+"\n"
    #print "FILTER SUMMARY: "+out
    outSegs = outBam.split("/")
    outDir = ""
    for i in range(len(outSegs)-1):
        outDir += outSegs[i] + "/"
    outFile = outDir+"filtering_result_summary.txt"
    #print "FILEDIR: "+outFile
    f = open(outFile, 'a')
    f.write(out)
    f.close()

def filter_sample_reads(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="filter_sample_reads v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _filter_sample_reads(logger, options)
