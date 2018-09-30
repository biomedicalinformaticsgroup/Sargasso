import os
import os.path
import schema
from . import constants


from . import file_writer as fw
from datetime import datetime

class MakefileWriter(object):
    def write_makefile(self, logger, options, sample_info):
        raise NotImplementedError()

class RnaseqMakefileWriter(MakefileWriter):
    pass

class ChipseqMakefileWriter(MakefileWriter):
    pass

class MakefileWriterManager(object):
    FileWriter = {"rnaseq" : RnaseqMakefileWriter,
            "chipseq": ChipseqMakefileWriter}

    def __init__(self, data_type, logger, options, sample_info):
        self.data_type = data_type
        self.logger = logger
        self.options = options
        self.sample_info = sample_info

    def get(self):
        mfw = self.FileWriter[self.data_type]
        return mfw(self.logger, self.options, self.sample_info)

