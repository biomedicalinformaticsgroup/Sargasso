import pysam


class SamUtils(object):

    def __init__(self):
        pass

    @staticmethod
    def open_samfile_for_read(filename):
        return pysam.Samfile(filename, "rb")

    @staticmethod
    def open_samfile_for_write(filename, template):
        return pysam.Samfile(filename, "wb", template=template)

    @staticmethod
    def all_hits(samfile):
        return samfile.fetch(until_eof=True)

    @staticmethod
    def hits_generator(samfile):
        last_hit_name = None
        current_hits = []

        for hit in SamUtils.all_hits(samfile):
            current_hit_name = hit.query_name

            if last_hit_name is None:
                last_hit_name = current_hit_name

            if current_hit_name == last_hit_name:
                current_hits.append(hit)
            else:
                yield current_hits
                last_hit_name = current_hit_name
                current_hits = [hit]

        yield current_hits

# class RnaseqSamUtils(SamUtils):
#     pass
#
#
# class ChipseqSamUtils(SamUtils):
#     pass
#
#
# class SamUtilsManager(Manager):
#     SAMUTILS = {'rnaseq': ChipseqSamUtils,
#                 'chipseq': RnaseqSamUtils}
#
#     @staticmethod
#     def get(data_type):
#         return SamUtilsManager.SAMUTILS[data_type](data_type)
