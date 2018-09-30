import docopt
import logging
import textwrap
import sys
import os.path
import schema

from schema import And, Or, Schema, Use
from datetime import datetime

from . import process
from . import file_writer as fw


from .__init__ import __version__
from . import filter_sample_reads

class Separator(object):

    DATA_TYPE=None

    DOC="""
    Usage:
        species_separator rnaseq [options]
        species_separator chipseq [options]
    Options:
        rnaseq run ss on rnaseq data
        chipseq run ss on chipseq data
    """

    def __init__(self, args):
        self.args = args

    def run(self):

        commandlineParser = CommandlineParserManager().get(self.DATA_TYPE)
        options = commandlineParser.parse(self.args, self.DOC)

        # Validate command-line options
        parameterValidator = ParameterValidatorManager().get(self.DATA_TYPE)
        sample_info = parameterValidator.validate(options)

        # Set up logger
        logger = LoggerManager(options).get()

        # Create output directory
        # os.mkdir(options[Option.OUTPUT_DIR])

        # Write Makefile to output directory
        makefileWriter = MakefileWriterManager().get(self.DATA_TYPE)
        makefileWriter.write(logger, options, sample_info)

        # Write Execution Record to output directory
        makefileWriter.write_execution_record(options)

        # If specified, execute the Makefile with nohup
        if options[Option.RUN_SEPARATION]:
            self._run_species_separation(logger, options)

        print(1111)
        return(0)

        # exec the make file
    def _run_species_separation(self,logger, options):
        """
        Executes the written Makefile with nohup.

        logger: logging object
        options: dictionary of command-line options
        """
        process.run_in_directory(options[Option.OUTPUT_DIR], "make")









class RnaseqSeparator(Separator):
    DOC="""
Usage:
    species_separator rnaseq 
    [--log-level=<log-level>]
    [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
    [--mismatch-threshold=<mismatch-threshold>]
    [--minmatch-threshold=<minmatch-threshold>]
    [--multimap-threshold=<multimap-threshold>]
    [--reject-multimaps]
    [--best] [--conservative] [--recall] [--permissive]
    [--run-separation]
    [--delete-intermediate]
    [--star-executable=<star-executable>]
    [--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>]
    <samples-file> <output-dir>
    (<species> <species-star-info>)
    (<species> <species-star-info>)
    ...

    
Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
    <samples-file>
        TSV file giving paths (relative to <reads-base-dir>) of raw RNA-seq read
        data files for each sample.
    <output-dir>
        Output directory into which Makefile will be written, and in which species
        separation will be performed.
    <species>
        Name of species.
    <species-star-info>
        Either a STAR index directory for the species, or a comma-separated list of
        (i) a GTF annotation file and (ii) a directory containing genome FASTA
        files for the species.
    <species-genome-fasta>
        Directory containing genome FASTA files for species.
    <species-index>
        STAR index directory for species.
    --reads-base-dir=<reads-base-dir>
        Base directory for raw RNA-seq read data files.
    -t <num-threads> --num-threads=<num-threads>
        Number of threads to use for parallel processing. [default: 1]
    --mismatch-threshold=<mismatch-threshold>
        Maximum percentage of bases allowed to be mismatches against the genome
        during filtering. For single-end reads, the total number of bases is the
        read length; for paired-end reads it is twice the read length [default: 0].
    --minmatch-threshold=<minmatch-threshold>
        Maximum percentage of read length allowed to not be mapped during
        filtering. Read length refers to the length of a single read in both
        single- and paired-end cases [default: 0].
    --multimap-threshold=<multimap-threshold>
        Maximum number of multiple mappings allowed during filtering [default: 1].
    --reject-multimaps
        If set, any read which multimaps to either species' genome will be rejected
        and not be assigned to either species.
    --best
        Adopt a filtering strategy that provides an excellent balance between
        sensitivity and specificity. Note that specifying this option overrides the
        values of the mismatch-threshold, minmatch-threshold and
        multimap-threshold options. In addition, reject-multimaps is turned off.
    --conservative
        Adopt a filtering strategy where minimising the number of reads
        mis-assigned to the wrong species takes foremost priority. Note that
        specifying this option overrides the values of the mismatch-threshold,
        minmatch-threshold and multimap-threshold options. In addition,
        reject-multimaps is turned on.
    --recall
        Adopt a filtering strategy where sensitivity is prioritised over
        specificity.  Note that specifying this option overrides the values of the
        mismatch-threshold, minmatch-threshold and multimap-threshold options. In
        addition, reject-multimaps is turned off.
    --permissive
        Adopt a filtering strategy where sensitivity is maximised. Note that
        specifying this option overrides the values of the mismatch-threshold,
        minmatch-threshold and multimap-threshold options. In addition,
        reject-multimaps is turned off.
    --run-separation
        If specified, species separation will be run; otherwise scripts to perform
        separation will be created but not run.
    --delete-intermediate
        Deletes the raw mapped BAMs and the sorted BAMs to free up space.
    --star-executable=<star-executable>
        Specify STAR executable path. Use this to run Sargasso with a particular
        version of STAR [default: STAR].
    --sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>
        Specify 'sambamba sort' temporary folder path [default: /tmp].
    
    Given a set of RNA-seq samples containing mixed-species read data, determine,
    where possible, from which of the species each read originated. Mapped
    reads are written to per-sample and -species specific output BAM files.
    
    Species separation of mixed-species RNA-seq data is performed in a number of
    stages, of which the most important steps are:
    
    1) Mapping of raw RNA-seq data to each species' genome using the STAR read
    aligner.
    2) Sorting of mapped RNA-seq data.
    3) Assignment of mapped, sorted RNA-seq reads to their correct species of
    origin.
    
    If the option "--run-separation" is not specified, a Makefile is written to the
    given output directory, via which all stages of species separation can be
    run under the user's control. If "--run-separation" is specified however, the
    Makefile is both written and executed, and all stages of species separation are
    performed automatically.
    
    n.b. Many stages of species separation can be executed across multiple threads
    by specifying the "--num-threads" option.
    
    e.g.:
    
    species_separator --reads-base-dir=/srv/data/rnaseq --num-threads 4 --run-separation samples.tsv my_results mouse /srv/data/genome/mouse/STAR_Index rat /srv/data/genome/rat/STAR_Index
    
"""
    DATA_TYPE = 'rnaseq'

class ChipseqSeparator(Separator):
    DOC="""not implemented"""
    DATA_TYPE = 'chipseq'


class CommandlineParser(object):
    def __init__(self):
        # self.DATA_TYPE=data_type
        # self.args = args
        # self.doc = doc
        # self.pv = ParameterValidator()
        pass

    def parse(self, args, doc):
        options = self._parse(args,doc)
        return options
        # self.pv.validate()

    def _parse(self,args,doc):
        # Read in command-line options
        docstring = self._substitute_common_options_into_usage(doc)
        options = docopt.docopt(docstring, argv=args, version="species_separator v" + __version__)
        return options

    def _substitute_common_options_into_usage(self, usage_msg, **substitutions):
        """
        Substitute common option and other interpolations into a usage message.

        Substitute help, version and logging level option specifications and
        descriptions into a script's usage message; also substitute any other
        interpolations specified by additional keyword arguments.
        usage_msg: A script's usage message.
        substitutions: Additional key=value interpolations to substitute into the
        usage message.
        """
        help_spec = "-h --help"
        help_desc = "Show this message."

        ver_spec = "-v --version"
        ver_desc = "Show version."

        log_spec = "{log_option}=<{log_level}>".format(
            log_option=Option._LOG_LEVEL_OPTION, log_level=Option.LOG_LEVEL)
        log_desc = ("Set logging level " +
                    "(one of {log_level_vals}) [default: info].").format(
            log_level_vals=str(Option.LEVELS.keys()))
        log_desc = '\n'.join(textwrap.wrap(
            log_desc, width=75, subsequent_indent="    "))

        return usage_msg.format(
            log_option_spec=log_spec, log_option_description=log_desc,
            help_option_spec=help_spec, help_option_description=help_desc,
            ver_option_spec=ver_spec, ver_option_description=ver_desc,
            **substitutions)



    def get_sample_info(self):
        raise NotImplementedError
        return sample_info

    def get_logger(self):
        raise NotImplementedError
        return self.logger



class RnaseqCommandlineParser(CommandlineParser):
    pass
class ChipseCommandlineParser(CommandlineParser):
    pass


class ParameterValidator(object):
    DATA_TYPE=None

    def __init__(self):
        pass

    def validate(self,options):
        raise NotImplementedError()

    def validate_dict_option(self, dict_option, values_dict, msg):
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

    def validate_dir_option(self, dir_option, msg, should_exist=True, nullable=False):
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
            validator = self._nullable_validator(validator)
        Schema(validator, error=msg).validate(dir_option)

    def validate_file_option(self, file_option, msg, should_exist=True, nullable=False):
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
            validator = self._nullable_validator(validator)
        Schema(validator, error=msg).validate(file_option)


    def validate_int_option(self, int_option, msg, min_val=None, nullable=False):
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
            validator = self._nullable_validator(validator)

        return Schema(validator, error=msg).validate(int_option)

    def validate_float_option(self, float_option, msg, min_val=None,
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
            validator = self._nullable_validator(validator)

        return Schema(validator, error=msg).validate(float_option)

    def check_boolean_value(self, option_string):
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


    def _nullable_validator(self, validator):
        return Or(validator, None)


class RnaseqParameterValidator(ParameterValidator):
    def validate(self,options):
        """
        Validate command line options are correctly specified.

        options: dictionary of command-line options.
        """
        try:
            self._validate_log_level(options)
            self.validate_dir_option(
                options[Option.READS_BASE_DIR], "Reads base directory does not exist",
                nullable=True)
            options[Option.NUM_THREADS] = self.validate_int_option(
                options[Option.NUM_THREADS],
                "Number of threads must be a positive integer",
                min_val=1, nullable=True)
            self.validate_file_option(
                options[Option.SAMPLES_FILE], "Could not open samples definition file")
            # debug remove comment in prod
            # opt.validate_dir_option(
            #     options[Option.OUTPUT_DIR], "Output directory should not exist",
            #     should_exist=False)


            if options[Option.OPTIMAL_STRATEGY]:
                options[Option.MISMATCH_THRESHOLD] = 1
                options[Option.MINMATCH_THRESHOLD] = 2
                options[Option.MULTIMAP_THRESHOLD] = 999999
                options[Option.REJECT_MULTIMAPS] = False
            elif options[Option.CONSERVATIVE_STRATEGY]:
                options[Option.MISMATCH_THRESHOLD] = 0
                options[Option.MINMATCH_THRESHOLD] = 0
                options[Option.MULTIMAP_THRESHOLD] = 1
                options[Option.REJECT_MULTIMAPS] = True
            elif options[Option.RECALL_STRATEGY]:
                options[Option.MISMATCH_THRESHOLD] = 2
                options[Option.MINMATCH_THRESHOLD] = 10
                options[Option.MULTIMAP_THRESHOLD] = 999999
                options[Option.REJECT_MULTIMAPS] = False
            elif options[Option.PERMISSIVE_STRATEGY]:
                options[Option.MISMATCH_THRESHOLD] = 25
                options[Option.MINMATCH_THRESHOLD] = 25
                options[Option.MULTIMAP_THRESHOLD] = 999999
                options[Option.REJECT_MULTIMAPS] = False

            # todo remove debug comment
            # filter_sample_reads.validate_threshold_options(
            #     options, Option.MISMATCH_THRESHOLD, Option.MINMATCH_THRESHOLD,
            #     Option.MULTIMAP_THRESHOLD)

            for i, species in enumerate(options[Option.SPECIES]):
                species_options = self._get_species_options(options, i)
                self._validate_species_options(species, species_options)

            sample_info = self._read_sample_info(options)
            # TODO: validate that all samples consistently have either single- or
            # paired-end reads
            self._validate_read_file(sample_info)

            return sample_info
        except schema.SchemaError as exc:
            exit("Exiting: " + exc.code)

    def _validate_read_file(self,sample_info):
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

    def _validate_log_level(self, options):
        """
        Check a command-line option specified logging level is valid.

        Check that the logging level specified through the command-line option
        log.LOG_LEVEL is one of the keys of the LEVELS dictionary.
        options: Dictionary mapping from command-line option strings to option
        values.
        """
        self.validate_dict_option(
            options[Option._LOG_LEVEL_OPTION], Option.LEVELS, "Invalid log level")

    def _read_sample_info(self, options):
        """
        Return an object encapsulating samples and their accompanying read files.
        options: dictionary of command-line options
        """
        sample_info = SampleInfo(options[Option.READS_BASE_DIR])

        for line in open(options[Option.SAMPLES_FILE], 'r'):
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

        species_options = { Option.SPECIES_NAME: options[Option.SPECIES][species_index] }

        star_infos = options[Option.SPECIES_STAR_INFO][species_index].split(",")

        if len(star_infos) == 1:
            species_options[Option.STAR_INDEX] = star_infos[0]
            species_options[Option.GTF_FILE] = None
            species_options[Option.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[Option.STAR_INDEX] = None
            species_options[Option.GTF_FILE] = star_infos[0]
            species_options[Option.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[Option.SPECIES_NAME]))

        return species_options

    def _validate_species_options(self, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        self.validate_file_option(
            species_options[Option.GTF_FILE],
            "Could not open species {species} GTF file".format(species=species),
            nullable=True)
        self.validate_dir_option(
            species_options[Option.GENOME_FASTA],
            "Genome FASTA directory for species {species} should exist".
                format(species=species),
            nullable=True)
        self.validate_dir_option(
            species_options[Option.STAR_INDEX],
            "STAR index directory for species {species} should exist".
                format(species=species),
            nullable=True)


class ChipseqParameterValidator(ParameterValidator):
    pass

class MakefileWriter(object):
    def __init__(self):
        pass

    def write(self):
        raise NotImplementedError()
    def write_execution_record(self):
        raise NotImplementedError()
class RnaseqMakefileWriter(MakefileWriter):

    # Write Makefile to output directory
    def write(self,logger, options, sample_info):
        print('Wring Rnaseq make file...')
        """
        Write Makefile to results directory to perform species separation.
        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        with fw.writing_to_file(
                fw.MakefileWriter, options[Option.OUTPUT_DIR], "Makefile") as writer:
            self._write_variable_definitions(logger, writer, options, sample_info)
            self._write_target_variable_definitions(logger, writer)
            self._write_phony_targets(logger, writer)
            self._write_all_target(logger, writer)
            self._write_filtered_reads_target(logger, writer, options)
            self._write_sorted_reads_target(logger, writer, options)
            self._write_mapped_reads_target(logger, writer, sample_info, options)
            # # self._write_masked_reads_target(logger, writer)
            self._write_collate_raw_reads_target(logger, writer, sample_info)
            # self._write_mask_star_index_targets(logger, writer, options)
            self._write_main_star_index_targets(logger, writer, options)
            self._write_clean_target(logger, writer)

    # Write Execution Record to output directory
    def write_execution_record(self,options):
        print('Wring Rnaseq execution_record...')
        """
        Write a log file containing all execution parameters in addition to the
        time and date of execution
        
        options: dictionary of command-line options
        """
        out_text = "Execution Record - {t}\n".format(
            t=str(datetime.now().isoformat()))
        out_text += "\n".join(["{desc}: {val}".format(
            desc=it[0], val=str(options[it[1]]))
            for it in Option.EXECUTION_RECORD_ENTRIES])

        out_file = os.path.join(options[Option.OUTPUT_DIR], "execution_record.txt")
        with open(out_file, 'w') as erf:
            erf.write(out_text)



    def _write_variable_definitions(self, logger, writer, options, sample_info):
        """
        Write variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        writer.set_variable(Option.NUM_THREADS_VARIABLE, options[Option.NUM_THREADS])
        writer.add_blank_line()

        writer.set_variable(Option.STAR_EXECUTABLE_VARIABLE, options[Option.STAR_EXECUTABLE])
        writer.add_blank_line()

        writer.set_variable(Option.SAMBAMBA_SORT_TMP_DIR_VARIABLE,options[Option.SAMBAMBA_SORT_TMP_DIR])
        writer.add_blank_line()

        sample_names = sample_info.get_sample_names()
        writer.set_variable(Option.SAMPLES_VARIABLE, " ".join(sample_names))
        writer.set_variable(
            Option.RAW_READS_DIRECTORY_VARIABLE,
            options[Option.READS_BASE_DIR] if options[Option.READS_BASE_DIR] else "/")
        writer.set_variable(
            Option.RAW_READS_LEFT_VARIABLE,
            " ".join([",".join(sample_info.get_left_reads(name))
                      for name in sample_names]))

        if sample_info.paired_end_reads():
            writer.set_variable(
                Option.RAW_READS_RIGHT_VARIABLE,
                " ".join([",".join(sample_info.get_right_reads(name))
                          for name in sample_names]))

        writer.add_blank_line()

        for i, species in enumerate(options[Option.SPECIES]):
            species_options = self._get_species_options(options, i)
            self._write_species_variable_definitions(
                logger, writer, species, species_options)

        writer.add_blank_line()

    def _write_target_variable_definitions(self, logger, writer):
        """
        Write target directory variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        writer.set_variable(Option.STAR_INDICES_TARGET, "star_indices")
        writer.set_variable(Option.COLLATE_RAW_READS_TARGET, "raw_reads")
        writer.set_variable(Option.MAPPED_READS_TARGET, "mapped_reads")
        writer.set_variable(Option.SORTED_READS_TARGET, "sorted_reads")
        writer.set_variable(Option.FILTERED_READS_TARGET, "filtered_reads")
        writer.add_blank_line()

    def _write_phony_targets(self, logger, writer):
        """
        Write phony target definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                ".PHONY", [Option.ALL_TARGET, Option.CLEAN_TARGET],
                raw_target=True, raw_dependencies=True):
            pass

    def _write_all_target(self, logger, writer):
        """
        Write main target definition to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                Option.ALL_TARGET, [Option.FILTERED_READS_TARGET], raw_target=True):
            pass

    def _write_filtered_reads_target(self, logger, writer, options):
        """
        Write target to separate reads by species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                Option.FILTERED_READS_TARGET, [Option.SORTED_READS_TARGET]):
            writer.add_comment(
                "For each sample, take the reads mapping to each genome and " +
                "filter them to their correct species of origin")
            writer.make_target_directory(Option.FILTERED_READS_TARGET)

            writer.add_command("filter_reads", [
                "\"{var}\"".format(var=writer.variable_val(Option.SAMPLES_VARIABLE)),
                writer.variable_val(Option.SORTED_READS_TARGET),
                writer.variable_val(Option.FILTERED_READS_TARGET),
                writer.variable_val(Option.NUM_THREADS_VARIABLE),
                options[Option.MISMATCH_THRESHOLD],
                options[Option.MINMATCH_THRESHOLD],
                options[Option.MULTIMAP_THRESHOLD],
                "--reject-multimaps" if options[Option.REJECT_MULTIMAPS] else "\"\"",
                "{sl}".format(sl=" ".join(options[Option.SPECIES]))])

            if options[Option.DELETE_INTERMEDIATE]:
                writer.remove_target_directory(Option.SORTED_READS_TARGET)


    def _write_sorted_reads_target(self, logger, writer, options):
        """
        Write target to sort reads by name to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        with writer.target_definition(
                Option.SORTED_READS_TARGET, [Option.MAPPED_READS_TARGET]):
            writer.add_comment(
                "For each sample, sort the mapped reads into read name order")
            writer.make_target_directory(Option.SORTED_READS_TARGET)

            writer.add_command(
                "sort_reads",
                ["\"{sl}\"".format(sl=" ".join(options[Option.SPECIES])),
                 "\"{var}\"".format(
                     var=writer.variable_val(Option.SAMPLES_VARIABLE)),
                 writer.variable_val(Option.NUM_THREADS_VARIABLE),
                 writer.variable_val(Option.MAPPED_READS_TARGET),
                 writer.variable_val(Option.SORTED_READS_TARGET),
                 writer.variable_val(Option.SAMBAMBA_SORT_TMP_DIR_VARIABLE)])

            if options[Option.DELETE_INTERMEDIATE]:
                writer.remove_target_directory(Option.MAPPED_READS_TARGET)

    def _write_mapped_reads_target(self, logger, writer, sample_info, options):
        """
        Write target to map reads to each species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """

        index_targets = ["{index}/{species}".format(
            index=writer.variable_val(Option.STAR_INDICES_TARGET),
            species=species)
            for species in options[Option.SPECIES]]

        index_targets.append(writer.variable_val(Option.COLLATE_RAW_READS_TARGET))

        with writer.target_definition(
                Option.MAPPED_READS_TARGET, index_targets, raw_dependencies=True):
            writer.add_comment(
                "Map reads for each sample to each species' genome")
            writer.make_target_directory(Option.MAPPED_READS_TARGET)

            map_reads_params = [
                "\"{sl}\"".format(sl=" ".join(options[Option.SPECIES])),
                "\"{var}\"".format(var=writer.variable_val(Option.SAMPLES_VARIABLE)),
                writer.variable_val(Option.STAR_INDICES_TARGET),
                writer.variable_val(Option.NUM_THREADS_VARIABLE),
                writer.variable_val(Option.COLLATE_RAW_READS_TARGET),
                writer.variable_val(Option.MAPPED_READS_TARGET)]

            map_reads_params.append(
                Option.PAIRED_END_READS_TYPE if sample_info.paired_end_reads()
                else Option.SINGLE_END_READS_TYPE)

            map_reads_params.append(options[Option.STAR_EXECUTABLE])

            writer.add_command("map_reads", map_reads_params)

    # def _write_masked_reads_target(logger, writer, options):
    # TODO: Write "masked_reads" target
    #pass

    def _write_collate_raw_reads_target(self, logger, writer, sample_info):
        """
        Write target to collect raw reads files to Makefile.

        logger: logging object
        writer: Makefile writer object
        sample_info: object encapsulating samples and their accompanying read files
        """
        with writer.target_definition(Option.COLLATE_RAW_READS_TARGET, []):
            writer.add_comment(
                "Create a directory with sub-directories for each sample, " +
                "each of which contains links to the input raw reads files " +
                "for that sample")
            writer.make_target_directory(Option.COLLATE_RAW_READS_TARGET)

            collate_raw_reads_params = [
                "\"{var}\"".format(var=writer.variable_val(Option.SAMPLES_VARIABLE)),
                writer.variable_val(Option.RAW_READS_DIRECTORY_VARIABLE),
                writer.variable_val(Option.COLLATE_RAW_READS_TARGET),
            ]

            if sample_info.paired_end_reads():
                collate_raw_reads_params += [
                    Option.PAIRED_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=writer.variable_val(Option.RAW_READS_LEFT_VARIABLE)),
                    "\"{var}\"".format(
                        var=writer.variable_val(Option.RAW_READS_RIGHT_VARIABLE))
                ]
            else:
                collate_raw_reads_params += [
                    Option.SINGLE_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=writer.variable_val(Option.RAW_READS_LEFT_VARIABLE)),
                    "\"\""
                ]

            writer.add_command("collate_raw_reads", collate_raw_reads_params)


    #def _write_mask_star_index_targets(logger, writer, options):
    # TODO: Write targets for STAR indices for mask sequences
    #pass


    def _write_species_main_star_index_target(self, logger, writer, species, species_options,star_executable):
        """
        Write target to create or link to STAR index for a species to Makefile.

        logger: logging object
        writer: Makefile writer object
        species_var: Makefile variable for species
        species_options: dictionary of command-line options for the particular
        species
        """
        target = "{index}/{spec}".format(
            index=writer.variable_val(Option.STAR_INDICES_TARGET), spec=species)

        with writer.target_definition(target, [], raw_target=True):
            if species_options[Option.GTF_FILE] is None:
                writer.make_target_directory(Option.STAR_INDICES_TARGET)
                writer.add_command(
                    "ln",
                    ["-s",
                     writer.variable_val(self._get_star_index_variable(species)),
                     target])
            else:
                writer.make_target_directory(target, raw_target=True)
                writer.add_command(
                    "build_star_index",
                    [writer.variable_val(self._get_genome_fasta_variable(species)),
                     writer.variable_val(self._get_gtf_file_variable(species)),
                     writer.variable_val(Option.NUM_THREADS_VARIABLE),
                     target,star_executable])

    def _write_main_star_index_targets(self, logger, writer, options):
        """
        Write targets to create or link to STAR indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        for i, species in enumerate(options[Option.SPECIES]):
            species_options = self._get_species_options(options, i)
            self._write_species_main_star_index_target(
                logger, writer, species, species_options,options[Option.STAR_EXECUTABLE])

    def _write_clean_target(self, logger, writer):
        """
        Write target to clean results directory to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(Option.CLEAN_TARGET, [], raw_target=True):
            writer.remove_target_directory(Option.STAR_INDICES_TARGET)
            writer.remove_target_directory(Option.COLLATE_RAW_READS_TARGET)
            writer.remove_target_directory(Option.MAPPED_READS_TARGET)
            writer.remove_target_directory(Option.SORTED_READS_TARGET)
            writer.remove_target_directory(Option.FILTERED_READS_TARGET)


    def _get_species_options(self, options, species_index):
        """
        Return a dictionary containing command-line options for a species.

        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = { Option.SPECIES_NAME: options[Option.SPECIES][species_index] }

        star_infos = options[Option.SPECIES_STAR_INFO][species_index].split(",")

        if len(star_infos) == 1:
            species_options[Option.STAR_INDEX] = star_infos[0]
            species_options[Option.GTF_FILE] = None
            species_options[Option.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[Option.STAR_INDEX] = None
            species_options[Option.GTF_FILE] = star_infos[0]
            species_options[Option.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[Option.SPECIES_NAME]))

        return species_options


    def _write_species_variable_definitions(self, logger, writer, species, species_options):
        """
        Write variable definitions for a particular species.

        logger: logging object
        writer: Makefile writer object
        species: species number string
        species_options: dictionary of options specific to a particular species.
        """
        if species_options[Option.GTF_FILE] is None:
            writer.set_variable(
                self._get_star_index_variable(species), species_options[Option.STAR_INDEX])
        else:
            writer.set_variable(
                self._get_gtf_file_variable(species), species_options[Option.GTF_FILE])
            writer.set_variable(
                self._get_genome_fasta_variable(species), species_options[Option.GENOME_FASTA])

    def _get_star_index_variable(self, species):
        """
        Return STAR index variable name for a species.

        species: species identifier
        """
        return "{species}_STAR_DIR".format(species=species.upper())


    def _get_gtf_file_variable(self, species):
        """
        Return GTF file variable name for a species.

        species: species identifier
        """
        return "{species}_GTF".format(species=species.upper())

    def _get_genome_fasta_variable(self, species):
        """
        Return genome FASTA variable name for a species.

        species: species identifier
        """
        return "{species}_GENOME_FASTA_DIR".format(species=species.upper())



class ChipseqMakefileWriter(MakefileWriter):
    pass


class Logger(object):

    def write_execution_record(self):
        raise NotImplementedError


class Manager(object):
    def _create(self):
        pass

    def get(self):
        pass

class LoggerManager(Manager):


    def __init__(self, options):
        self.options = options
        self.logger = self._create(options)

    def _create(self,options):
        """
        Return a Logger instance with a command-line specified severity threshold.

        Return a Logger instance with severity threshold specified by the command
        line option log.LOG_LEVEL. Log messages will be written to standard out.
        options: Dictionary mapping from command-line option strings to option
        values.
        """
        return self._get_logger(sys.stderr, options[Option._LOG_LEVEL_OPTION])

    def _get_logger(self, stream, level):
        """
        Return a Logger instance with the specified severity threshold.

        Return a Logger instance with the specified severity threshold, where the
        threshold level should be a key of the 'LEVELS' dictionary. Log messages
        will contain the current time and message severity level.
        stream: Output stream to which the logger will write messages.
        level: Severity threshold level, which should be a key of the 'LEVELS'
        dictionary.
        """
        formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s: %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')
        handler = logging.StreamHandler(stream)
        handler.setFormatter(formatter)
        logger = logging.getLogger(__name__)
        logger.setLevel(Option.LEVELS[level])
        logger.addHandler(handler)
        return logger

    def get(self):
        return self.logger



class SeparatorManager(object):
    SEPARATORS = {"rnaseq": RnaseqSeparator,
                  "chipseq": ChipseqSeparator}

    def __init__(self,data_type,args):
        self.args = args
        self.data_type = data_type
        self.separator = self._creat()

    def _creat(self):
        separator = self.SEPARATORS[self.data_type]
        return separator(self.args)

    def get(self):
        return self.separator


class CommandlineParserManager(Manager):
    PARSERS = {"rnaseq" : RnaseqCommandlineParser,
               "chipseq": ChipseCommandlineParser}

    def __init__(self):
        # self.args = args
        # self.data_type = data_type
        # self.doc = doc
        # self.parser=self._create()
        pass

    def _create(self,data_type):
        commandlineParser = self.PARSERS[data_type]
        return commandlineParser()

    def get(self,data_type):
        return self._create(data_type)

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

class MakefileWriterManager(Manager):
    mkw = {"rnaseq" : RnaseqMakefileWriter,
           "chipseq": ChipseqMakefileWriter}

    def _create(self,data_type):
        return self.mkw[data_type]()

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



class Option(object):
    SAMPLES_FILE = "<samples-file>"
    OUTPUT_DIR = "<output-dir>"
    SPECIES = "<species>"
    SPECIES_STAR_INFO = "<species-star-info>"
    READS_BASE_DIR = "--reads-base-dir"
    NUM_THREADS = "--num-threads"
    MISMATCH_THRESHOLD = "--mismatch-threshold"
    MINMATCH_THRESHOLD = "--minmatch-threshold"
    MULTIMAP_THRESHOLD = "--multimap-threshold"
    REJECT_MULTIMAPS = "--reject-multimaps"
    OPTIMAL_STRATEGY = "--best"
    CONSERVATIVE_STRATEGY = "--conservative"
    RECALL_STRATEGY = "--recall"
    PERMISSIVE_STRATEGY = "--permissive"
    RUN_SEPARATION = "--run-separation"
    DELETE_INTERMEDIATE = "--delete-intermediate"
    STAR_EXECUTABLE = "--star-executable"
    SAMBAMBA_SORT_TMP_DIR = "--sambamba-sort-tmp-dir"

    SPECIES_NAME = "species-name"
    GTF_FILE = "gtf-file"
    GENOME_FASTA = "genome-fasta"
    STAR_INDEX = "star-index"

    ALL_TARGET = "all"
    CLEAN_TARGET = "clean"
    STAR_INDICES_TARGET = "STAR_INDICES"
    COLLATE_RAW_READS_TARGET = "COLLATE_RAW_READS"
    MAPPED_READS_TARGET = "MAPPED_READS"
    SORTED_READS_TARGET = "SORTED_READS"
    FILTERED_READS_TARGET = "FILTERED_READS"

    NUM_THREADS_VARIABLE = "NUM_THREADS"
    STAR_EXECUTABLE_VARIABLE = "STAR_EXECLEVELSUTABLE"
    SAMBAMBA_SORT_TMP_DIR_VARIABLE = "SAMBAMBA_SORT_TMP_DIR"
    SAMPLES_VARIABLE = "SAMPLES"
    RAW_READS_DIRECTORY_VARIABLE = "RAW_READS_DIRECTORY"
    RAW_READS_LEFT_VARIABLE = "RAW_READS_FILES_1"
    RAW_READS_RIGHT_VARIABLE = "RAW_READS_FILES_2"

    SINGLE_END_READS_TYPE = "single"
    PAIRED_END_READS_TYPE = "paired"

    LOG_LEVEL = "log-level"
    LEVELS = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL
    }
    _LOG_LEVEL_OPTION = "--" + LOG_LEVEL

    EXECUTION_RECORD_ENTRIES = [
        ["Samples File", SAMPLES_FILE],
        ["Output Dir", OUTPUT_DIR],
        ["Species", SPECIES],
        ["Species STAR info", SPECIES_STAR_INFO],
        ["Reads Base Dir", READS_BASE_DIR],
        ["Number of Threads", NUM_THREADS],
        ["Mismatch Threshold", MISMATCH_THRESHOLD],
        ["Minmatch Threshold", MINMATCH_THRESHOLD],
        ["Multimap Threshold", MULTIMAP_THRESHOLD],
        ["Reject Multimaps", REJECT_MULTIMAPS],
        ["Optimal Strategy", OPTIMAL_STRATEGY],
        ["Conservative Strategy", CONSERVATIVE_STRATEGY],
        ["Recall Strategy", RECALL_STRATEGY],
        ["Run Separation", RUN_SEPARATION],
        ["Delete Intermediate" , DELETE_INTERMEDIATE],
        ["Star Executable Path", STAR_EXECUTABLE],
        ["Sambamba Sort Tmp Dir", SAMBAMBA_SORT_TMP_DIR],
    ]




