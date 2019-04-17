import schema

import sargasso.filter.filter_controllers as fcs
import sargasso.filter.sample_filterer as sfs
import sargasso.separator.file_writer as fw
import sargasso.separator.commandline_parser as clp
import sargasso.separator.parameter_validator as pv
import sargasso.separator.separators as sps


_DATA_TYPE_OPTION = "<data-type>"
_DATA_TYPES = {}


def _data_type(cls):
    _DATA_TYPES[cls.NAME] = cls()
    return cls


def get_data_type_manager(args, docstring):
    data_type = clp.CommandlineParser.parse_datatype(
        args, docstring, _DATA_TYPE_OPTION)

    try:
        pv.ParameterValidator.validate_dict_option(
            data_type, _DATA_TYPES, "Invalid data type.")
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)

    return _DATA_TYPES[data_type]


class _BaseDataTypeManager(object):
    def __init__(
            self,
            command_line_parser_cls,
            parameter_validator_cls,
            makefile_writer_cls,
            separator_cls,
            filter_controller_cls,
            sample_filterer_cls):

        self.command_line_parser = command_line_parser_cls()
        self.parameter_validator = parameter_validator_cls()
        self.makefile_writer = makefile_writer_cls(self.NAME)
        self.separator = separator_cls(
                self.command_line_parser,
                self.parameter_validator,
                self.makefile_writer,
                fw.ExecutionRecordWriter())
        self.filter_controller = filter_controller_cls(
                self.NAME,
                self.command_line_parser)
        self.sample_filterer = sample_filterer_cls(
                self.command_line_parser)

    def get_command_line_parser(self):
        return self.command_line_parser

    def get_parameter_validator(self):
        return self.parameter_validator

    def get_makefile_writer(self):
        return self.makefile_writer

    def get_separator(self):
        return self.separator

    def get_filter_controller(self):
        return self.filter_controller

    def get_sample_filterer(self):
        return self.sample_filterer


@_data_type
class _RnaSeqDataTypeManager(_BaseDataTypeManager):
    NAME = "rnaseq"

    def __init__(self):
        _BaseDataTypeManager.__init__(
                self,
                clp.RnaSeqCommandlineParser,
                pv.RnaSeqParameterValidator,
                fw.RnaSeqMakefileWriter,
                sps.RnaSeqSeparator,
                fcs.RnaSeqFilterController,
                sfs.RnaSeqSampleFilterer)


@_data_type
class _DnaSeqDataTypeManager(_BaseDataTypeManager):
    NAME = "dnaseq"

    def __init__(self):
        _BaseDataTypeManager.__init__(
                self,
                clp.DnaSeqCommandlineParser,
                pv.DnaSeqParameterValidator,
                fw.DnaSeqMakefileWriter,
                sps.DnaSeqSeparator,
                fcs.DnaSeqFilterController,
                sfs.DnaSeqSampleFilterer)


@_data_type
class _BisulfiteDataTypeManager(_BaseDataTypeManager):
    NAME = "bisulfite"

    def __init__(self):
        _BaseDataTypeManager.__init__(
            self,
            clp.BisulfiteCommandlineParser,
            pv.BisulfiteParameterValidator,
            fw.BisulfiteMakefileWriter,
            sps.BisulfiteSeparator,
            fcs.BisulfiteFilterController,
            sfs.BisulfiteSampleFilterer)

