import pysam


def open_samfile_for_read(filename):
    return pysam.Samfile(filename, "rb")


def open_samfile_for_write(filename, template):
    return pysam.Samfile(filename, "wb", template=template)


def all_hits(samfile):
    return samfile.fetch(until_eof=True)


def get_read_length(hit):
    return hit.query_length


def get_total_length(primary_hits):
    total_length = 0
    for hit in primary_hits:
        total_length += hit.query_length
    return total_length


def get_multimaps(hit):
    return hit.get_tag("NH")


def get_mismatches(hit):
    return hit.get_tag("nM")


def get_alignment_scores(hit):
    return hit.get_tag("AS")


def get_cigar(hit):
    return hit.cigartuples

def is_primary(hit):
    return not hit.is_secondary


def hits_generator(samfile):
    last_hit_name = None
    current_hits = None

    for hit in all_hits(samfile):
        if current_hits is None:
            current_hits = []
            
        current_hit_name = hit.query_name

        if last_hit_name is None:
            last_hit_name = current_hit_name

        if current_hit_name == last_hit_name:
            current_hits.append(hit)
        else:
            yield current_hits
            last_hit_name = current_hit_name
            current_hits = [hit]

        if current_hits is not None:
            yield current_hits
