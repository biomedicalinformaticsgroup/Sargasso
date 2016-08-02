import pysam


def open_samfile_for_read(filename):
    return pysam.Samfile(filename, "rb")


def open_samfile_for_write(filename, template):
    return pysam.Samfile(filename, "wb", template=template)


def all_hits(samfile):
    return samfile.fetch(until_eof=True)


def get_read_length(hit):
    return hit.query_length


def get_total_length(hit):
    total_length = hit.query_length
    if hit.is_paired:
        total_length *= 2
    return total_length


def get_multimaps(hit):
    return hit.get_tag("NH")


def get_mismatches(hit):
    return hit.get_tag("nM")


def get_alignment_scores(hit):
    return hit.get_tag("AS")


def get_cigar(hit):
    return hit.cigarstring


def hits_generator(samfile):
    last_hit_name = None
    current_hits = []

    for hit in all_hits(samfile):

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
