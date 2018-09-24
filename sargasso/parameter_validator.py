from . import options
from . import options as opt
from .__init__ import __version__
import schema
import os
import docopt
from . import constants
from . import filter_sample_reads

class Validator(object):

    def __init__(self,args):
        self.args = args

    def parse_command_line_options(self,args,doc):
        # Read in command-line options
        docstring = opt.substitute_common_options_into_usage(doc)
        options = docopt.docopt(docstring, argv=args, version="species_separator v" + __version__)
        return options

    def validate_command_line_options(self,options):
         raise NotImplementedError()



class RnaseqValidator(Validator):

    def validate_command_line_options(self,options):
        """
        Validate command line options are correctly specified.

        options: dictionary of command-line options.
        """

        try:
            opt.validate_log_level(options)
            opt.validate_dir_option(
                options[constants.READS_BASE_DIR], "Reads base directory does not exist",
                nullable=True)
            options[constants.NUM_THREADS] = opt.validate_int_option(
                options[constants.NUM_THREADS],
                "Number of threads must be a positive integer",
                min_val=1, nullable=True)
            opt.validate_file_option(
                options[constants.SAMPLES_FILE], "Could not open samples definition file")
            # debug remove comment in prod
            # opt.validate_dir_option(
            #     options[constants.OUTPUT_DIR], "Output directory should not exist",
            #     should_exist=False)


            if options[constants.OPTIMAL_STRATEGY]:
                options[constants.MISMATCH_THRESHOLD] = 1
                options[constants.MINMATCH_THRESHOLD] = 2
                options[constants.MULTIMAP_THRESHOLD] = 999999
                options[constants.REJECT_MULTIMAPS] = False
            elif options[constants.CONSERVATIVE_STRATEGY]:
                options[constants.MISMATCH_THRESHOLD] = 0
                options[constants.MINMATCH_THRESHOLD] = 0
                options[constants.MULTIMAP_THRESHOLD] = 1
                options[constants.REJECT_MULTIMAPS] = True
            elif options[constants.RECALL_STRATEGY]:
                options[constants.MISMATCH_THRESHOLD] = 2
                options[constants.MINMATCH_THRESHOLD] = 10
                options[constants.MULTIMAP_THRESHOLD] = 999999
                options[constants.REJECT_MULTIMAPS] = False
            elif options[constants.PERMISSIVE_STRATEGY]:
                options[constants.MISMATCH_THRESHOLD] = 25
                options[constants.MINMATCH_THRESHOLD] = 25
                options[constants.MULTIMAP_THRESHOLD] = 999999
                options[constants.REJECT_MULTIMAPS] = False

            filter_sample_reads.validate_threshold_options(
                options, constants.MISMATCH_THRESHOLD, constants.MINMATCH_THRESHOLD,
                constants.MULTIMAP_THRESHOLD)

            for i, species in enumerate(options[constants.SPECIES]):
                species_options = self._get_species_options(options, i)
                self._validate_species_options(species, species_options)

            sample_info = self._read_sample_info(options)
            # TODO: validate that all samples consistently have either single- or
            # paired-end reads
            sample_info.validate()

            return sample_info
        except schema.SchemaError as exc:
            exit("Exiting: " + exc.code)


    def _get_species_options(self, options, species_index):
        """
        Return a dictionary containing command-line options for a species.
        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = { constants.SPECIES_NAME: options[constants.SPECIES][species_index] }

        star_infos = options[constants.SPECIES_STAR_INFO][species_index].split(",")

        if len(star_infos) == 1:
            species_options[constants.STAR_INDEX] = star_infos[0]
            species_options[constants.GTF_FILE] = None
            species_options[constants.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[constants.STAR_INDEX] = None
            species_options[constants.GTF_FILE] = star_infos[0]
            species_options[constants.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[constants.SPECIES_NAME]))

        return species_options

    def _validate_species_options(self, species, species_options):
        """
        Validate command-line options for a species are correctly specified.
        species: species identification string
        species_options: dictionary of options specific to a particular species.
        """
        opt.validate_file_option(
            species_options[constants.GTF_FILE],
            "Could not open species {species} GTF file".format(species=species),
            nullable=True)
        opt.validate_dir_option(
            species_options[constants.GENOME_FASTA],
            "Genome FASTA directory for species {species} should exist".
                format(species=species),
            nullable=True)
        opt.validate_dir_option(
            species_options[constants.STAR_INDEX],
            "STAR index directory for species {species} should exist".
                format(species=species),
            nullable=True)

    def _read_sample_info(self, options):
        """
        Return an object encapsulating samples and their accompanying read files.
        options: dictionary of command-line options
        """
        sample_info = SampleInfo(options[constants.READS_BASE_DIR])

        for line in open(options[constants.SAMPLES_FILE], 'r'):
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


class ChipseqValidator(Validator):
    def validate_command_line_options(self,options):
        raise NotImplementedError()


class ValidatorManager(object):
    VALIDATORS = {"rnaseq" : RnaseqValidator,
                  "chipseq": ChipseqValidator}

    def __init__(self, data_type, args):
        self.args = args
        self.data_type = data_type
        self.validator = self.creat_validator(data_type,args)

    def creat_validator(self,data_type,args):
        validator = self.VALIDATORS[data_type]
        return validator(args)

    def get_validator(self):
        return self.validator



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

    def validate(self):
        """
        Validate all raw reads files exist.
        """
        for reads_set in [self.left_reads, self.right_reads]:
            for reads_file_list in reads_set.values():
                for reads_file in reads_file_list:
                    self.validate_read_file(reads_file)

    def validate_read_file(self, reads_file):
        """
        Validate a particular raw reads file exists.
        reads_file: Path to a particular raw reads file (possibly truncated,
        and needing to be joined with 'base_reads_dir').
        """
        if self.base_reads_dir:
            reads_file = os.path.join(self.base_reads_dir, reads_file)
        opt.validate_file_option(
            reads_file, "Could not open reads file")