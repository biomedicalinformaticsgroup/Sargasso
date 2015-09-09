import pysam


def open_samfile_for_read(filename):
    return pysam.Samfile(filename, "rb")


def open_samfile_for_write(filename, template):
    return pysam.Samfile(filename, "wb", template=template)


def all_reads(samfile):
    return samfile.fetch(until_eof=True)


def get_edit_distance(read):
    return read.opt("NM")


def get_alignment_score(read):
    return read.opt("AS")
