import contextlib
import textwrap
import os
import os.path
import sargasso.separator.options as opts

from sargasso.utils import log
from datetime import datetime


class Writer(object):
    def __init__(self):
        self.lines = []

    def _add_line(self, line_string):
        self.lines.append(line_string)

    @contextlib.contextmanager
    def writing_to_file(self, directory, filename):
        try:
            yield
        finally:
            self._write_to_file(directory, filename)

    def _write_to_file(self, directory, filename):
        with open(os.path.join(directory, filename), "w") as output_file:
            output_file.write("\n".join(self.lines) + '\n')

    def write(self, *args):
        raise NotImplementedError('Need to implement in subclass')


class MakefileWriter(Writer):
    ALL_TARGET = "all"
    CLEAN_TARGET = "clean"
    MAPPER_INDICES_TARGET = "MAPPER_INDICES"
    COLLATE_RAW_READS_TARGET = "COLLATE_RAW_READS"
    MAPPED_READS_TARGET = "MAPPED_READS"
    SORTED_READS_TARGET = "SORTED_READS"
    FILTERED_READS_TARGET = "FILTERED_READS"

    DATA_TYPE_VARIABLE = "DATA_TYPE"
    NUM_THREADS_VARIABLE = "NUM_THREADS"
    MAPPER_EXECUTABLE_VARIABLE = "MAPPER_EXECUTABLE"
    SAMBAMBA_SORT_TMP_DIR_VARIABLE = "SAMBAMBA_SORT_TMP_DIR"
    SAMPLES_VARIABLE = "SAMPLES"
    RAW_READS_DIRECTORY_VARIABLE = "RAW_READS_DIRECTORY"
    RAW_READS_LEFT_VARIABLE = "RAW_READS_FILES_1"
    RAW_READS_RIGHT_VARIABLE = "RAW_READS_FILES_2"

    SINGLE_END_READS_TYPE = "single"
    PAIRED_END_READS_TYPE = "paired"

    INDENT = '\t'

    def __init__(self, data_type):
        Writer.__init__(self)
        self.data_type = data_type
        self.indent_level = 0

    def indent(self):
        self.indent_level += 1

    def deindent(self):
        self.indent_level -= 1

    def add_line(self, line_string):
        line_string = MakefileWriter.INDENT * self.indent_level + line_string
        Writer._add_line(self, line_string)

    def add_blank_line(self):
        self.add_line("")

    @contextlib.contextmanager
    def target_definition(self, target, dependencies,
                          raw_target=False, raw_dependencies=False):

        self.add_line("{tar}: {deps}".format(
            tar=self.variable_val(target, raw_target),
            deps=" ".join([self.variable_val(d, raw_dependencies)
                for d in dependencies])))
        self.indent()

        try:
            yield
        finally:
            self.deindent()
            self.add_blank_line()

    def add_comment(self, comment):
        lines = textwrap.wrap(
            comment, initial_indent="# ", subsequent_indent="# ",
            width=75 - self.indent_level * len(MakefileWriter.INDENT))
        for line in lines:
            self.add_line(line)

    def set_variable(self, variable, value):
        self.add_line("{var}={val}".format(var=variable, val=value))

    @classmethod
    def variable_val(cls, variable, raw=False):
        return variable if raw else "$({var})".format(var=variable)

    def make_target_directory(self, target, raw_target=False):
        self.add_command("mkdir", ["-p", self.variable_val(target, raw_target)])

    def remove_target_directory(self, target, raw_target=False):
        self.add_command("rm", ["-rf", self.variable_val(target, raw_target)])

    def add_command(self, command_name, options):
        line_elements = [command_name] + options
        self.add_line(" ".join([str(l) for l in line_elements]))

    def _write_variable_definitions(self, options, sample_info, species_options):
        """
        Write variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        self.set_variable(MakefileWriter.NUM_THREADS_VARIABLE,
                          options[opts.NUM_THREADS])
        self.add_blank_line()

        self.set_variable(MakefileWriter.DATA_TYPE_VARIABLE,
                          options[opts.DATA_TYPE_ARG])
        self.add_blank_line()

        self.set_variable(MakefileWriter.MAPPER_EXECUTABLE_VARIABLE,
                          options[opts.MAPPER_EXECUTABLE])
        self.add_blank_line()

        self.set_variable(MakefileWriter.SAMBAMBA_SORT_TMP_DIR_VARIABLE,
                          options[opts.SAMBAMBA_SORT_TMP_DIR])
        self.add_blank_line()

        sample_names = sample_info.get_sample_names()
        self.set_variable(MakefileWriter.SAMPLES_VARIABLE, " ".join(sample_names))
        self.set_variable(
            MakefileWriter.RAW_READS_DIRECTORY_VARIABLE,
            options[opts.READS_BASE_DIR] if options[opts.READS_BASE_DIR] else "/")
        self.set_variable(
            MakefileWriter.RAW_READS_LEFT_VARIABLE,
            " ".join([",".join(sample_info.get_left_reads(name))
                      for name in sample_names]))

        if sample_info.paired_end_reads():
            self.set_variable(
                MakefileWriter.RAW_READS_RIGHT_VARIABLE,
                " ".join([",".join(sample_info.get_right_reads(name))
                          for name in sample_names]))

        self.add_blank_line()

        for i, species in enumerate(options[opts.SPECIES_ARG]):
            self._write_species_variable_definitions(species, species_options[i])

        self.add_blank_line()

    def _write_species_variable_definitions(self, species, species_options):
        raise NotImplementedError()

    def _write_target_variable_definitions(self):
        """
        Write target directory variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        self.set_variable(MakefileWriter.MAPPER_INDICES_TARGET, "mapper_indexes")
        self.set_variable(MakefileWriter.COLLATE_RAW_READS_TARGET, "raw_reads")
        self.set_variable(MakefileWriter.MAPPED_READS_TARGET, "mapped_reads")
        self.set_variable(MakefileWriter.SORTED_READS_TARGET, "sorted_reads")
        self.set_variable(MakefileWriter.FILTERED_READS_TARGET, "filtered_reads")
        self.add_blank_line()

    def _write_phony_targets(self):
        """
        Write phony target definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.target_definition(
                ".PHONY", [MakefileWriter.ALL_TARGET, MakefileWriter.CLEAN_TARGET],
                raw_target=True, raw_dependencies=True):
            pass

    def _write_all_target(self):
        """
        Write main target definition to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.target_definition(
                MakefileWriter.ALL_TARGET,
                [MakefileWriter.FILTERED_READS_TARGET],
                raw_target=True):
            pass

    def _write_filtered_reads_target(self, options):
        """
        Write target to separate reads by species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.target_definition(
                MakefileWriter.FILTERED_READS_TARGET,
                [MakefileWriter.SORTED_READS_TARGET]):
            self.add_comment(
                "For each sample, take the reads mapping to each genome and " +
                "filter them to their correct species of origin")
            self.make_target_directory(MakefileWriter.FILTERED_READS_TARGET)

            self.add_command("filter_reads", [
                self.variable_val(MakefileWriter.DATA_TYPE_VARIABLE),
                "\"{var}\"".format(
                    var=self.variable_val(MakefileWriter.SAMPLES_VARIABLE)),
                self.variable_val(MakefileWriter.SORTED_READS_TARGET),
                self.variable_val(MakefileWriter.FILTERED_READS_TARGET),
                self.variable_val(MakefileWriter.NUM_THREADS_VARIABLE),
                options[opts.MISMATCH_THRESHOLD],
                options[opts.MINMATCH_THRESHOLD],
                options[opts.MULTIMAP_THRESHOLD],
                "--reject-multimaps" if options[opts.REJECT_MULTIMAPS] else "\"\"",
                options[log.LOG_LEVEL_OPTION],
                "{sl}".format(sl=" ".join(options[opts.SPECIES_ARG]))])

            if options[opts.DELETE_INTERMEDIATE]:
                self.remove_target_directory(MakefileWriter.SORTED_READS_TARGET)

    def _write_sorted_reads_target(self, options):
        """
        Write target to sort reads by name to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        with self.target_definition(
                MakefileWriter.SORTED_READS_TARGET,
                [MakefileWriter.MAPPED_READS_TARGET]):
            self.add_comment(
                "For each sample, sort the mapped reads into read name order")
            self.make_target_directory(MakefileWriter.SORTED_READS_TARGET)

            self.add_command(
                "sort_reads",
                ["\"{sl}\"".format(sl=" ".join(options[opts.SPECIES_ARG])),
                 "\"{var}\"".format(
                     var=self.variable_val(MakefileWriter.SAMPLES_VARIABLE)),
                 self.variable_val(MakefileWriter.NUM_THREADS_VARIABLE),
                 self.variable_val(MakefileWriter.MAPPED_READS_TARGET),
                 self.variable_val(MakefileWriter.SORTED_READS_TARGET),
                 self.variable_val(MakefileWriter.SAMBAMBA_SORT_TMP_DIR_VARIABLE)])

            if options[opts.DELETE_INTERMEDIATE]:
                self.remove_target_directory(MakefileWriter.MAPPED_READS_TARGET)

    def _write_mapped_reads_target(self, sample_info, options):
        """
        Write target to map reads to each species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """

        index_targets = ["{index}/{species}".format(
            index=self.variable_val(MakefileWriter.MAPPER_INDICES_TARGET),
            species=species)
            for species in options[opts.SPECIES_ARG]]

        index_targets.append(
                self.variable_val(MakefileWriter.COLLATE_RAW_READS_TARGET))

        with self.target_definition(
                MakefileWriter.MAPPED_READS_TARGET, index_targets,
                raw_dependencies=True):

            self.add_comment(
                "Map reads for each sample to each species' genome")
            self.make_target_directory(MakefileWriter.MAPPED_READS_TARGET)

            map_reads_params = \
                ["\"{sl}\"".format(sl=" ".join(options[opts.SPECIES_ARG])),
                 "\"{var}\"".format(
                     var=self.variable_val(MakefileWriter.SAMPLES_VARIABLE)),
                 self.variable_val(MakefileWriter.MAPPER_INDICES_TARGET),
                 self.variable_val(MakefileWriter.NUM_THREADS_VARIABLE),
                 self.variable_val(MakefileWriter.COLLATE_RAW_READS_TARGET),
                 self.variable_val(MakefileWriter.MAPPED_READS_TARGET),
                 MakefileWriter.PAIRED_END_READS_TYPE if
                     sample_info.paired_end_reads() else
                     MakefileWriter.SINGLE_END_READS_TYPE,
                 options[opts.MAPPER_EXECUTABLE]]

            self.add_command("map_reads_" + self.data_type, map_reads_params)

    def _write_collate_raw_reads_target(self, sample_info):
        """
        Write target to collect raw reads files to Makefile.

        logger: logging object
        writer: Makefile writer object
        sample_info: object encapsulating samples and their accompanying read files
        """
        with self.target_definition(MakefileWriter.COLLATE_RAW_READS_TARGET, []):
            self.add_comment(
                "Create a directory with sub-directories for each sample, " +
                "each of which contains links to the input raw reads files " +
                "for that sample")
            self.make_target_directory(MakefileWriter.COLLATE_RAW_READS_TARGET)

            collate_raw_reads_params = [
                "\"{var}\"".format(
                    var=self.variable_val(MakefileWriter.SAMPLES_VARIABLE)),
                self.variable_val(MakefileWriter.RAW_READS_DIRECTORY_VARIABLE),
                self.variable_val(MakefileWriter.COLLATE_RAW_READS_TARGET),
            ]

            if sample_info.paired_end_reads():
                collate_raw_reads_params += [
                    MakefileWriter.PAIRED_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=self.variable_val(
                            MakefileWriter.RAW_READS_LEFT_VARIABLE)),
                    "\"{var}\"".format(
                        var=self.variable_val(
                            MakefileWriter.RAW_READS_RIGHT_VARIABLE))
                ]
            else:
                collate_raw_reads_params += [
                    MakefileWriter.SINGLE_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=self.variable_val(
                            MakefileWriter.RAW_READS_LEFT_VARIABLE)),
                    "\"\""
                ]

            self.add_command("collate_raw_reads", collate_raw_reads_params)

    @classmethod
    def _get_species_options(cls, options, species_index):
        return options['species_options'][species_index]

    def write(self, *args):
        raise NotImplementedError()


class RnaSeqMakefileWriter(MakefileWriter):

    # Write Makefile to output directory
    def write(self, logger, options):
        """
        Write Makefile to results directory to perform species separation.
        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        sample_info = options[opts.SAMPLE_INFO_INDEX]
        species_options = options[opts.SPECIES_OPTIONS_INDEX]

        with self.writing_to_file(options[opts.OUTPUT_DIR_ARG], "Makefile"):
            self._write_variable_definitions(options, sample_info, species_options)
            self._write_target_variable_definitions()
            self._write_phony_targets()
            self._write_all_target()
            self._write_filtered_reads_target(options)
            self._write_sorted_reads_target(options)
            self._write_mapped_reads_target(sample_info, options)
            # # self._write_masked_reads_target(logger)
            self._write_collate_raw_reads_target(sample_info)
            # self._write_mask_star_index_targets(logger, options)
            self._write_main_star_index_targets(options)
            self._write_clean_target()

    # def _write_masked_reads_target(logger, options):
    # TODO: Write "masked_reads" target
    # pass

    # def _write_mask_star_index_targets(logger, options):
    # TODO: Write targets for STAR indices for mask sequences
    # pass

    def _write_main_star_index_targets(self, options):
        """
        Write targets to create or link to STAR indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        for i, species in enumerate(options[opts.SPECIES_ARG]):
            species_options = self._get_species_options(options, i)
            self._write_species_main_star_index_target(
                species, species_options, options[opts.MAPPER_EXECUTABLE])

    def _write_species_main_star_index_target(
            self, species, species_options, executable):
        """
        Write target to create or link to STAR index for a species to Makefile.

        logger: logging object
        writer: Makefile writer object
        species_var: Makefile variable for species
        species_options: dictionary of command-line options for the particular
        species
        """
        target = "{index}/{spec}".format(
            index=self.variable_val(MakefileWriter.MAPPER_INDICES_TARGET),
            spec=species)

        with self.target_definition(target, [], raw_target=True):
            if species_options[opts.GTF_FILE] is None:
                self.make_target_directory(MakefileWriter.MAPPER_INDICES_TARGET)
                self.add_command(
                    "ln",
                    ["-s",
                     self.variable_val(self._get_star_index_variable(species)),
                     target])
            else:
                self.make_target_directory(target, raw_target=True)
                self.add_command(
                    "build_star_index",
                    [self.variable_val(self._get_genome_fasta_variable(species)),
                     self.variable_val(self._get_gtf_file_variable(species)),
                     self.variable_val(MakefileWriter.NUM_THREADS_VARIABLE),
                     target, executable])

    def _write_clean_target(self):
        """
        Write target to clean results directory to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.target_definition(MakefileWriter.CLEAN_TARGET, [],
                                    raw_target=True):
            self.remove_target_directory(MakefileWriter.MAPPER_INDICES_TARGET)
            self.remove_target_directory(MakefileWriter.COLLATE_RAW_READS_TARGET)
            self.remove_target_directory(MakefileWriter.MAPPED_READS_TARGET)
            self.remove_target_directory(MakefileWriter.SORTED_READS_TARGET)
            self.remove_target_directory(MakefileWriter.FILTERED_READS_TARGET)

    def _write_species_variable_definitions(self, species, species_options):
        """
        Write variable definitions for a particular species.

        logger: logging object
        writer: Makefile writer object
        species: species number string
        species_options: dictionary of options specific to a particular species.
        """
        if species_options[opts.GTF_FILE] is None:
            self.set_variable(
                self._get_star_index_variable(species),
                species_options[opts.MAPPER_INDEX])
        else:
            self.set_variable(
                self._get_gtf_file_variable(species),
                species_options[opts.GTF_FILE])
            self.set_variable(
                self._get_genome_fasta_variable(species),
                species_options[opts.GENOME_FASTA])

    @classmethod
    def _get_star_index_variable(cls, species):
        """
        Return STAR index variable name for a species.

        species: species identifier
        """
        return "{species}_STAR_DIR".format(species=species.upper())

    @classmethod
    def _get_gtf_file_variable(cls, species):
        """
        Return GTF file variable name for a species.

        species: species identifier
        """
        return "{species}_GTF".format(species=species.upper())

    @classmethod
    def _get_genome_fasta_variable(cls, species):
        """
        Return genome FASTA variable name for a species.

        species: species identifier
        """
        return "{species}_GENOME_FASTA_DIR".format(species=species.upper())


class DnaSeqMakefileWriter(MakefileWriter):
    # Write Makefile to output directory
    def write(self, logger, options):
        """
        Write Makefile to results directory to perform species separation.
        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        # todo ask Owen this
        sample_info = options[opts.SAMPLE_INFO_INDEX]
        species_options = options[opts.SPECIES_OPTIONS_INDEX]
        with self.writing_to_file(options[opts.OUTPUT_DIR_ARG], "Makefile"):
            self._write_variable_definitions(options, sample_info, species_options)
            self._write_target_variable_definitions()
            self._write_phony_targets()
            self._write_all_target()
            self._write_filtered_reads_target(options)
            self._write_sorted_reads_target(options)
            self._write_mapped_reads_target(sample_info, options)
            # # # self._write_masked_reads_target(logger)
            self._write_collate_raw_reads_target(sample_info)
            # # self._write_mask_star_index_targets(logger, options)
            self._write_main_bowtie2_index_targets(options)
            # self._write_clean_target(logger)

    def _write_main_bowtie2_index_targets(self, options):
        """
        Write targets to create or link to Bowtie indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        for i, species in enumerate(options[opts.SPECIES_ARG]):
            species_options = self._get_species_options(options, i)
            self._write_species_main_bowtie2_index_target(
                species, species_options, options[opts.MAPPER_INDEX_EXECUTABLE])

    def _write_species_main_bowtie2_index_target(
            self, species, species_options, executable):
        """
        Write target to create or link to Bowtie index for a species to Makefile.

        logger: logging object
        writer: Makefile writer object
        species_var: Makefile variable for species
        species_options: dictionary of command-line options for the particular
        species
        """
        target = "{index}/{spec}".format(
            index=self.variable_val(MakefileWriter.MAPPER_INDICES_TARGET),
            spec=species)

        with self.target_definition(target, [], raw_target=True):
            if species_options[opts.GENOME_FASTA] is None:
                self.make_target_directory(MakefileWriter.MAPPER_INDICES_TARGET)
                self.add_command(
                    "ln",
                    ["-s",
                     self.variable_val(self._get_bowtie2_index_variable(species)),
                     target])
            else:
                self.make_target_directory(target, raw_target=True)
                self.add_command(
                    "build_bowtie2_index",
                    [self.variable_val(self._get_genome_fasta_variable(species)),
                     self.variable_val(MakefileWriter.NUM_THREADS_VARIABLE),
                     target, executable])

    def _write_species_variable_definitions(self, species, species_options):
        """
        Write variable definitions for a particular species.

        logger: logging object
        writer: Makefile writer object
        species: species number string
        species_options: dictionary of options specific to a particular species.
        """
        if species_options[opts.GENOME_FASTA] is None:
            self.set_variable(
                self._get_bowtie2_index_variable(species),
                species_options[opts.MAPPER_INDEX])
        else:
            self.set_variable(
                self._get_genome_fasta_variable(species),
                species_options[opts.GENOME_FASTA])

    @classmethod
    def _get_bowtie2_index_variable(cls, species):
        return "{species}_BOWTIE2_DIR".format(species=species.upper())

    @classmethod
    def _get_genome_fasta_variable(cls, species):
        return "{species}_GENOME_FASTA_FILE".format(species=species.upper())


class ExecutionRecordWriter(Writer):
    EXECUTION_RECORD_ENTRIES = [
        ["Data Type", opts.DATA_TYPE_ARG],
        ["Samples File", opts.SAMPLES_FILE_ARG],
        ["Output Dir", opts.OUTPUT_DIR_ARG],
        ["Species", opts.SPECIES_ARG],
        ["Species info", opts.SPECIES_INFO_ARG],
        ["Reads Base Dir", opts.READS_BASE_DIR],
        ["Number of Threads", opts.NUM_THREADS],
        ["Mismatch Threshold", opts.MISMATCH_THRESHOLD],
        ["Minmatch Threshold", opts.MINMATCH_THRESHOLD],
        ["Multimap Threshold", opts.MULTIMAP_THRESHOLD],
        ["Reject Multimaps", opts.REJECT_MULTIMAPS],
        ["Optimal Strategy", opts.OPTIMAL_STRATEGY],
        ["Conservative Strategy", opts.CONSERVATIVE_STRATEGY],
        ["Recall Strategy", opts.RECALL_STRATEGY],
        ["Permissive Strategy", opts.PERMISSIVE_STRATEGY],
        ["Run Separation", opts.RUN_SEPARATION],
        ["Delete Intermediate", opts.DELETE_INTERMEDIATE],
        ["Mapper Executable Path", opts.MAPPER_EXECUTABLE],
        ["Mapper Indexing Executable Path", opts.MAPPER_INDEX_EXECUTABLE],
        ["Sambamba Sort Tmp Dir", opts.SAMBAMBA_SORT_TMP_DIR],
    ]

    def write(self, options):
        """
        Write a log file containing all execution parameters in addition to the
        time and date of execution

        options: dictionary of command-line options
        """
        out_text = "Execution Record - {t}\n".format(
            t=str(datetime.now().isoformat()))
        out_text += "\n".join(["{desc}: {val}".format(
            desc=it[0], val=str(options[it[1]]))
            for it in self.EXECUTION_RECORD_ENTRIES])

        out_file = os.path.join(options[opts.OUTPUT_DIR_ARG],
                                "execution_record.txt")

        with open(out_file, 'w') as erf:
            erf.write(out_text)
