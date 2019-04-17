import os.path
import schema
import sargasso.separator.options as opts

from schema import And, Or, Schema, Use
from sargasso.utils import log


class ParameterValidator(object):

    @classmethod
    def validate(cls, options):
        """
        Validate command line options are correctly specified.

        options: dictionary of command-line options.
        """
        sample_info = options[opts.SAMPLE_INFO_INDEX]
        species_options = options[opts.SPECIES_OPTIONS_INDEX]

        try:
            cls.validate_log_level(options)
            cls.validate_dir_option(
                options[opts.READS_BASE_DIR],
                "Reads base directory does not exist",
                nullable=True)
            options[opts.NUM_THREADS] = cls.validate_int_option(
                options[opts.NUM_THREADS],
                "Number of threads must be a positive integer",
                min_val=1, nullable=True)
            cls.validate_file_option(
                options[opts.SAMPLES_FILE_ARG],
                "Could not open samples definition file")
            cls.validate_dir_option(
                options[opts.OUTPUT_DIR_ARG],
                "Output directory should not exist",
                should_exist=False)

            cls.validate_threshold_options(
                options, opts.MISMATCH_THRESHOLD, opts.MINMATCH_THRESHOLD,
                opts.MULTIMAP_THRESHOLD)

            for i, species in enumerate(options[opts.SPECIES_ARG]):
                cls._validate_species_options(species, species_options[i])

            # TODO: validate that all samples consistently have either single- or
            # paired-end reads
            cls._validate_read_file(sample_info)

            return sample_info
        except schema.SchemaError as exc:
            exit("Exiting: " + exc.code)

    @classmethod
    def validate_log_level(cls, options):
        """
        Check a command-line option specified logging level is valid.

        Check that the logging level specified through the command-line option
        log.LOG_LEVEL is one of the keys of the LEVELS dictionary.
        options: Dictionary mapping from command-line option strings to option
        values.
        """
        ParameterValidator.validate_dict_option(
            options[log.LOG_LEVEL_OPTION], log.LEVELS, "Invalid log level")

    @classmethod
    def validate_dict_option(cls, dict_option, values_dict, msg):
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

    @classmethod
    def validate_dir_option(cls, dir_option, msg, should_exist=True, nullable=False):
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

    @classmethod
    def validate_file_option(cls, file_option, msg,
                             should_exist=True, nullable=False):
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

    @classmethod
    def validate_int_option(cls, int_option, msg, min_val=None, nullable=False):
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

    @classmethod
    def validate_float_option(cls, float_option, msg, min_val=None,
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

    @classmethod
    def check_boolean_value(cls, option_string):
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

    @classmethod
    def _nullable_validator(cls, validator):
        return Or(validator, None)

    @classmethod
    def _validate_read_file(cls, sample_info):
        """
        Validate all raw reads files exist.
        """
        for reads_set in [sample_info.left_reads, sample_info.right_reads]:
            for reads_file_list in reads_set.values():
                for reads_file in reads_file_list:
                    if sample_info.base_reads_dir:
                        reads_file = os.path.join(sample_info.base_reads_dir,
                                                  reads_file)
                    cls.validate_file_option(
                        reads_file, "Could not open reads file")

    @classmethod
    def _validate_species_options(cls, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        raise NotImplementedError()

    @classmethod
    def validate_threshold_options(
        cls, options, mismatch_opt_name, minmatch_opt_name, multimap_opt_name):

        options[mismatch_opt_name] = ParameterValidator.validate_float_option(
            options[mismatch_opt_name],
            "Maximum percentage of mismatches must be a float between 0 and 100",
            0, 100, True)
        options[minmatch_opt_name] = ParameterValidator.validate_float_option(
            options[minmatch_opt_name],
            "Maximum percentage of read length which does not match must be a " +
            "float between 0 and 100", 0, 100, True)
        options[multimap_opt_name] = ParameterValidator.validate_int_option(
            options[multimap_opt_name],
            "Maximum number of multiple mappings must be a positive integer",
            1, True)


class RnaSeqParameterValidator(ParameterValidator):

    @classmethod
    def _validate_species_options(cls, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        cls.validate_file_option(
            species_options[opts.GTF_FILE],
            "Could not open species {species} GTF file".format(species=species),
            nullable=True)
        cls.validate_dir_option(
            species_options[opts.GENOME_FASTA],
            "Genome FASTA directory for species {species} should exist".format(
                species=species),
            nullable=True)
        cls.validate_dir_option(
            species_options[opts.MAPPER_INDEX],
            "STAR index directory for species {species} should exist".format(
                species=species),
            nullable=True)


class DnaSeqParameterValidator(ParameterValidator):

    @classmethod
    def _validate_species_options(cls, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        cls.validate_file_option(
            species_options[opts.GENOME_FASTA],
            "Could not open species {species} FASTA file".format(species=species),
            nullable=True)
        cls.validate_dir_option(
            species_options[opts.MAPPER_INDEX],
            "Bowtie2 index directory for species {species} should exist".format(
                species=species),
            nullable=True)

class BisulfiteParameterValidator(DnaSeqParameterValidator):
    @classmethod
    def _validate_species_options(cls, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        cls.validate_file_option(
            species_options[opts.GENOME_FASTA],
            "Could not open species {species} FASTA file".format(species=species),
            nullable=True)
        cls.validate_dir_option(
            species_options[opts.MAPPER_INDEX],
            "Bisulfite index directory for species {species} should exist".format(
                species=species),
            nullable=True)
