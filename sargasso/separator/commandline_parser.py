import textwrap
import docopt
import schema
import sargasso.separator.options as opts

from sargasso import __version__
from sargasso.utils import log


class CommandlineParser(object):
    @classmethod
    def parse(cls, args, doc):
        docstring = cls._substitute_common_options_into_usage(doc)
        return docopt.docopt(docstring, argv=args,
                             version="species_separator v" + __version__,
                             options_first=False)

    @classmethod
    def parse_parameters(cls, args, doc):
        options = cls.parse(args, doc)
        options[opts.SAMPLE_INFO_INDEX] = cls.parse_sample_info(options)
        options[opts.SPECIES_OPTIONS_INDEX] = cls.parse_species_options(options)
        options = cls._parse_sargasso_strategy(options)
        return options

    @classmethod
    def parse_datatype(cls, args, doc, data_type_string):
        docstring = CommandlineParser._substitute_common_options_into_usage(doc)
        options = docopt.docopt(docstring, argv=args,
                                version="species_separator v" + __version__,
                                options_first=True)
        return options[data_type_string]

    @classmethod
    def _parse_sargasso_strategy(cls, options):
        if options[opts.OPTIMAL_STRATEGY]:
            options[opts.MISMATCH_THRESHOLD] = 1
            options[opts.MINMATCH_THRESHOLD] = 2
            options[opts.MULTIMAP_THRESHOLD] = 999999
            options[opts.REJECT_MULTIMAPS] = False
        elif options[opts.CONSERVATIVE_STRATEGY]:
            options[opts.MISMATCH_THRESHOLD] = 0
            options[opts.MINMATCH_THRESHOLD] = 0
            options[opts.MULTIMAP_THRESHOLD] = 1
            options[opts.REJECT_MULTIMAPS] = True
        elif options[opts.RECALL_STRATEGY]:
            options[opts.MISMATCH_THRESHOLD] = 2
            options[opts.MINMATCH_THRESHOLD] = 10
            options[opts.MULTIMAP_THRESHOLD] = 999999
            options[opts.REJECT_MULTIMAPS] = False
        elif options[opts.PERMISSIVE_STRATEGY]:
            options[opts.MISMATCH_THRESHOLD] = 25
            options[opts.MINMATCH_THRESHOLD] = 25
            options[opts.MULTIMAP_THRESHOLD] = 999999
            options[opts.REJECT_MULTIMAPS] = False
        return options

    @classmethod
    def _substitute_common_options_into_usage(cls, usage_msg, **substitutions):
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
            log_option=log.LOG_LEVEL_OPTION, log_level=log.LOG_LEVEL)
        log_desc = ("Set logging level " +
                    "(one of {log_level_vals}) [default: info].").format(
            log_level_vals=str(log.LEVELS.keys()))
        log_desc = '\n'.join(textwrap.wrap(
            log_desc, width=75, subsequent_indent="    "))

        return usage_msg.format(
            log_option_spec=log_spec, log_option_description=log_desc,
            help_option_spec=help_spec, help_option_description=help_desc,
            ver_option_spec=ver_spec, ver_option_description=ver_desc,
            **substitutions)

    @classmethod
    def parse_sample_info(cls, options):
        """
        Return an object encapsulating samples and their accompanying read files.
        options: dictionary of command-line options
        """
        sample_info = SampleInfo(options[opts.READS_BASE_DIR])

        for line in open(options[opts.SAMPLES_FILE_ARG], 'r'):
            sample_data = line.split()

            if len(sample_data) < 2 or len(sample_data) > 3:
                line = line.rstrip('\n')
                if len(line) > 80:
                    line = line[0:77] + "..."
                raise schema.SchemaError(
                    None, "Sample file line should contain sample name and " +
                          "lists of read files (or first and second pairs of " +
                          "read files), separated by whitespace: " +
                          "\n{info}".format(info=line))

            sample_info.add_sample_data(sample_data)

        return sample_info

    @classmethod
    def parse_species_options(cls, options):
        species_options = {}
        for i, species in enumerate(options[opts.SPECIES_ARG]):
            species_options[i] = cls._parse_species_options(options, i)
        return species_options

    @classmethod
    def _parse_species_options(cls, options, species_index):
        raise NotImplementedError('Need to implement in subclass')


class RnaSeqCommandlineParser(CommandlineParser):

    @classmethod
    def _parse_species_options(cls, options, species_index):
        """
        Return a dictionary containing command-line options for a species.
        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = {opts.SPECIES_NAME:
                           options[opts.SPECIES_ARG][species_index]}

        star_infos = options[opts.SPECIES_INFO_ARG][species_index].split(",")

        if len(star_infos) == 1:
            species_options[opts.MAPPER_INDEX] = star_infos[0]
            species_options[opts.GTF_FILE] = None
            species_options[opts.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[opts.MAPPER_INDEX] = None
            species_options[opts.GTF_FILE] = star_infos[0]
            species_options[opts.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[opts.SPECIES_NAME]))

        return species_options


class DnaSeqCommandlineParser(CommandlineParser):
    @classmethod
    def _parse_species_options(cls, options, species_index):
        """
        Return a dictionary containing command-line options for a species.
        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = {opts.SPECIES_NAME:
                           options[opts.SPECIES_ARG][species_index]}

        bowtie_infos = options[opts.SPECIES_INFO_ARG][species_index]

        if options[opts.SPECIES_INFO_ARG][species_index].endswith('.fa'):
            species_options[opts.MAPPER_INDEX] = None
            species_options[opts.GENOME_FASTA] = bowtie_infos
        else:
            species_options[opts.MAPPER_INDEX] = bowtie_infos
            species_options[opts.GENOME_FASTA] = None

        return species_options


class BisulfiteCommandlineParser(DnaSeqCommandlineParser):
    #todo bisulfite
    pass


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
