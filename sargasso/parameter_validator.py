from schema import And, Or, Schema, Use
from options import Options
import os.path
from factory import Manager
import schema




class ParameterValidator(object):
    DATA_TYPE=None

    def __init__(self):
        pass

    def validate(self, options):
        raise NotImplementedError()

    @staticmethod
    def validate_log_level(options):
        """
        Check a command-line option specified logging level is valid.

        Check that the logging level specified through the command-line option
        log.LOG_LEVEL is one of the keys of the LEVELS dictionary.
        options: Dictionary mapping from command-line option strings to option
        values.
        """
        ParameterValidator.validate_dict_option(
            options[Options.LOG_LEVEL_OPTION], Options.LEVELS, "Invalid log level")


    @staticmethod
    def validate_dict_option( dict_option, values_dict, msg):
        """
        Check if a command line option is a dictionary key.

        Check if a command line option is a key in the specified dictionary and, if
        so, return the corresponding dictionary value. If the option is not a key
        in the dictionary, a SchemaError is raised.

        dict_option: The command line option, a string.
        values_dict: A dictionary, in which a valid command line option should be a
        key.
        msg: Text for the SchemaError exception raised if the test fails.
        """
        msg = "{msg}: '{opt}'.".format(msg=msg, opt=dict_option)
        return Schema(Use(lambda x: values_dict[x]), error=msg). \
            validate(dict_option)
    @staticmethod
    def validate_dir_option(dir_option, msg, should_exist=True, nullable=False):
        """
        Check if a directory specified by a command line option exists or not.

        Check if the directory specified by a command line option exists or not, as
        indicated by the parameter 'should_exist'. The option can be allowed to
        equal 'None' if 'nullable' is set to True. If the specified directory
        fails the test, a SchemaError is raised.

        dir_option: A string, the path to the directory.
        msg: Text for the SchemaError exception raised if the test fails.
        should_exist: Determines if the directory is checked for existence of
        non-existence.
        nullable: If set to True, the command line option is allowed to be 'None'
        (i.e. the option has not been specified).
        """
        msg = "{msg}: '{dir}'.".format(msg=msg, dir=dir_option)
        validator = os.path.isdir if should_exist else \
            lambda f: not os.path.exists(f)
        if nullable:
            validator = ParameterValidator._nullable_validator(validator)
        Schema(validator, error=msg).validate(dir_option)

    @staticmethod
    def validate_file_option( file_option, msg, should_exist=True, nullable=False):
        """
        Check if a file specified by a command line option exists or not.

        Check if the file specified by a command line option exists or not, as
        indicated by the parameter 'should_exist'. The option can be allowed to
        equal 'None' if 'nullable' is set to True. If the specified file fails the
        test, a SchemaError is raised.

        file_option: A string, the path to the file.
        msg: Text for the SchemaError exception raised if the test fails.
        should_exist: Determines if the file is checked for existence or
        non-existence.
        """
        msg = "{msg}: '{file}'.".format(msg=msg, file=file_option)
        validator = open if should_exist else \
            lambda f: not os.path.exists(f)
        if nullable:
            validator = ParameterValidator._nullable_validator(validator)
        Schema(validator, error=msg).validate(file_option)

    @staticmethod
    def validate_int_option(int_option, msg, min_val=None, nullable=False):
        """
        Check if a command line option is an integer.

        Check if a command line option string represents a valid integer and, if
        so, return the integer value. If 'min_val' is specified, the integer must
        be greater than or equal to the given minimum value. The option can be
        allowed to equal 'None' if 'nullable' is set to True. If the option is not
        a valid integer, a SchemaError is raised.

        int_option: The command line option, a string.
        msg: Text for the SchemaError exception raised if the test fails.
        min_val: If set, the integer must be greater than or equal to this value.
        nullable: If set to True, the command line option is allowed to be 'None'
        (i.e. the option has not been specified).
        """
        msg = "{msg}: '{val}'".format(msg=msg, val=int_option)
        validator = Use(int)
        if min_val is not None:
            validator = And(validator, lambda x: x >= min_val)
        if nullable:
            validator = ParameterValidator._nullable_validator(validator)

        return Schema(validator, error=msg).validate(int_option)

    @staticmethod
    def validate_float_option( float_option, msg, min_val=None,
                              max_val=None, nullable=False):
        """
        Check if a command line option is a floating point number.

        Check if a command line option string represents a valid floating point
        number and, if so, return the float value. If 'min_val' is specified, the
        float must be greater than or equal to the given minimum value. If the
        option is not a valid float, a SchemaError is raised.

        float_option: The command line option, a string.
        msg: Text for the SchemaError exception raised if the test fails.
        min_val: If set, the float must be greater than or equal to this value.
        max_val: If set, the float must be less than or equal to this value.
        nullable: If set to True, the command line option is allowed to be 'None'
        (i.e. the option has not been specified).
        """
        msg = "{msg}: '{val}'".format(msg=msg, val=float_option)
        validator = Use(float)
        if min_val is not None:
            validator = And(validator, lambda x: x >= min_val)
        if max_val is not None:
            validator = And(validator, lambda x: x <= max_val)
        if nullable:
            validator = ParameterValidator._nullable_validator(validator)

        return Schema(validator, error=msg).validate(float_option)

    @staticmethod
    def check_boolean_value( option_string):
        """
        Validates that a command line option string represents a boolean value.

        Check if a command line option string represents a valid boolean value.
        Accepted strings for True are "true", "t", "yes", "y" or any cased variants
        thereof. Accepted string for False are "false", "f", "no", "n" or any cased
        variants thereof. Any other values are considered invalid, and supplying
        them will cause an exception to be raised.

        option_string: A command line option string representing a boolean value.
        """
        option_string = option_string.lower()
        if option_string in ["true", "t", "yes", "y"]:
            return True
        elif option_string in ["false", "f", "no", "n"]:
            return False
        else:
            raise Exception("Can't convert '{o}' to bool.".format(o=option_string))

    @staticmethod
    def _nullable_validator(validator):
        return Or(validator, None)

    @staticmethod
    def validate_datatype(data_type):
        ParameterValidator.validate_dict_option(
            data_type, Options.SUPPORTED_DATATYPE, "Invalid data type.")


class RnaseqParameterValidator(ParameterValidator):
    def validate(self,options):
        """
        Validate command line options are correctly specified.

        options: dictionary of command-line options.
        """
        try:
            self.validate_log_level(options)
            self.validate_dir_option(
                options[Options.READS_BASE_DIR], "Reads base directory does not exist",
                nullable=True)
            options[Options.NUM_THREADS] = self.validate_int_option(
                options[Options.NUM_THREADS],
                "Number of threads must be a positive integer",
                min_val=1, nullable=True)
            self.validate_file_option(
                options[Options.SAMPLES_FILE], "Could not open samples definition file")
            # debug remove comment in prod
            # opt.validate_dir_option(
            #     options[Options.OUTPUT_DIR], "Output directory should not exist",
            #     should_exist=False)


            if options[Options.OPTIMAL_STRATEGY]:
                options[Options.MISMATCH_THRESHOLD] = 1
                options[Options.MINMATCH_THRESHOLD] = 2
                options[Options.MULTIMAP_THRESHOLD] = 999999
                options[Options.REJECT_MULTIMAPS] = False
            elif options[Options.CONSERVATIVE_STRATEGY]:
                options[Options.MISMATCH_THRESHOLD] = 0
                options[Options.MINMATCH_THRESHOLD] = 0
                options[Options.MULTIMAP_THRESHOLD] = 1
                options[Options.REJECT_MULTIMAPS] = True
            elif options[Options.RECALL_STRATEGY]:
                options[Options.MISMATCH_THRESHOLD] = 2
                options[Options.MINMATCH_THRESHOLD] = 10
                options[Options.MULTIMAP_THRESHOLD] = 999999
                options[Options.REJECT_MULTIMAPS] = False
            elif options[Options.PERMISSIVE_STRATEGY]:
                options[Options.MISMATCH_THRESHOLD] = 25
                options[Options.MINMATCH_THRESHOLD] = 25
                options[Options.MULTIMAP_THRESHOLD] = 999999
                options[Options.REJECT_MULTIMAPS] = False

            # todo remove debug comment
            # filter_sample_reads.validate_threshold_options(
            #     options, Options.MISMATCH_THRESHOLD, Options.MINMATCH_THRESHOLD,
            #     Options.MULTIMAP_THRESHOLD)

            for i, species in enumerate(options[Options.SPECIES]):
                species_options = self._get_species_options(options, i)
                self._validate_species_options(species, species_options)

            sample_info = self._read_sample_info(options)
            # TODO: validate that all samples consistently have either single- or
            # paired-end reads
            self._validate_read_file(sample_info)

            return sample_info
        except schema.SchemaError as exc:
            exit("Exiting: " + exc.code)

    def _validate_read_file(self, sample_info):
        """
        Validate all raw reads files exist.
        """
        for reads_set in [sample_info.left_reads, sample_info.right_reads]:
            for reads_file_list in reads_set.values():
                for reads_file in reads_file_list:
                    if sample_info.base_reads_dir:
                        reads_file = os.path.join(sample_info.base_reads_dir, reads_file)
                    self.validate_file_option(
                        reads_file, "Could not open reads file")

    def _read_sample_info(self, options):
        """
        Return an object encapsulating samples and their accompanying read files.
        options: dictionary of command-line options
        """
        sample_info = SampleInfo(options[Options.READS_BASE_DIR])

        for line in open(options[Options.SAMPLES_FILE], 'r'):
            sample_data = line.split()

            if len(sample_data) < 2 or len(sample_data) > 3:
                line = line.rstrip('\n')
                if len(line) > 80:
                    line = line[0:77] + "..."
                raise schema.SchemaError(
                    None, "Sample file line should contain sample name and " +
                          "lists of read files (or first and second pairs of read " +
                          "files), separated by whitespace: \n{info}".format(info=line))

            sample_info.add_sample_data(sample_data)

        return sample_info

    def _get_species_options(self, options, species_index):
        """
        Return a dictionary containing command-line options for a species.
        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = { Options.SPECIES_NAME: options[Options.SPECIES][species_index] }

        star_infos = options[Options.SPECIES_STAR_INFO][species_index].split(",")

        if len(star_infos) == 1:
            species_options[Options.STAR_INDEX] = star_infos[0]
            species_options[Options.GTF_FILE] = None
            species_options[Options.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[Options.STAR_INDEX] = None
            species_options[Options.GTF_FILE] = star_infos[0]
            species_options[Options.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[Options.SPECIES_NAME]))

        return species_options

    def _validate_species_options(self, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        self.validate_file_option(
            species_options[Options.GTF_FILE],
            "Could not open species {species} GTF file".format(species=species),
            nullable=True)
        self.validate_dir_option(
            species_options[Options.GENOME_FASTA],
            "Genome FASTA directory for species {species} should exist".
                format(species=species),
            nullable=True)
        self.validate_dir_option(
            species_options[Options.STAR_INDEX],
            "STAR index directory for species {species} should exist".
                format(species=species),
            nullable=True)

class ChipseqParameterValidator(ParameterValidator):
    pass

class ParameterValidatorManager(Manager):
    VALIDATORS = {"rnaseq" : RnaseqParameterValidator,
                  "chipseq": ChipseqParameterValidator}

    def __init__(self):
        # self.data_type = data_type
        # self.options = options
        # self.validator = self._create()
        pass

    def _create(self,data_type):
        parameterValidator = self.VALIDATORS[data_type]
        return parameterValidator()

    def get(self,data_type):
        return self._create(data_type)


class SampleInfo(object):

    """
    Encapsulates sample names and their accompanying reads files.
    """

    def __init__(self, base_reads_dir):
        """
        Create object.
        base_reads_dir: Common base path to all raw reads files.
        """
        self.base_reads_dir = base_reads_dir
        self.left_reads = {}
        self.right_reads = {}
        self.paired_end = None

    def get_sample_names(self):
        """
        Return a list of sample names.
        """
        return self.left_reads.keys()

    def get_left_reads(self, sample_name):
        """
        Return list of first in pair reads files for sample.
        """
        return self.left_reads[sample_name]

    def get_right_reads(self, sample_name):
        """
        Return list of second in pair reads files for sample.
        """
        return self.right_reads[sample_name]

    def paired_end_reads(self):
        """
        Returns True iff sample read files are paired-end data.
        """
        return self.paired_end

    def add_sample_data(self, sample_data):
        """
        Add information for a single sample.
        sample_data: list of strings - sample name, comma-separated list of
        read files (or comma-separated list of first in pair reads files,
        comma-separated list of second in pair reads files).
        """
        sample_name = sample_data[0]
        self.left_reads[sample_name] = sample_data[1].split(',')

        if len(sample_data) == 3:
            self.right_reads[sample_name] = sample_data[2].split(',')
            self.paired_end = True
        else:
            self.paired_end = False






