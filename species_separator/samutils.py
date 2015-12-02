import pysam


def open_samfile_for_read(filename):
    return pysam.Samfile(filename, "rb")


def open_samfile_for_write(filename, template):
    return pysam.Samfile(filename, "wb", template=template)


def all_hits(samfile):
    return samfile.fetch(until_eof=True)


def get_multimaps(hit):
    return hit.get_tag("NH")


def get_mismatches(hit):
    return hit.get_tag("nM")

def get_alignment_scores(hit):
    return hit.get_tag("AS")

def get_cigar(hit):
    return hit.cigarstring
