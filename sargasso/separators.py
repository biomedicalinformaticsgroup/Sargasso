import docopt
from . import options as opt
import os
from . import process
import schema
from . import constants
from .__init__ import __version__
from . import parameter_validator as pv
from . import file_writer as fw
from datetime import datetime




class Separator(object):
    DOC="""
    Usage:
        species_separator rnaseq [options]
        species_separator chipseq [options]
    Options:
        rnaseq run ss on rnaseq data
        chipseq run ss on chipseq data
    """

    def __init__(self,data_type, args):
        self.args = args
        self.data_type = data_type
        self.validator = pv.ValidatorManager(self.data_type, args).get_validator()


     # Write Makefile to output directory
    def _write_makefile(self,logger, options, sample_info):
        raise NotImplementedError()

    # Write Execution Record to output directory
    def _write_execution_record(self,options):
        raise NotImplementedError()

    # exec the make file
    def _run_species_separation(self,logger, options):
        """
        Executes the written Makefile with nohup.

        logger: logging object
        options: dictionary of command-line options
        """
        process.run_in_directory(options[constants.OUTPUT_DIR], "make")

    def run(self):

        options = self.validator.parse_command_line_options(self.args,self.DOC)

        # Validate command-line options
        sample_info = self.validator.validate_command_line_options(options)

        # Set up logger
        logger = opt.get_logger_for_options(options)

        # Create output directory
        #os.mkdir(options[constants.OUTPUT_DIR])

        # Write Makefile to output directory
        self._write_makefile(logger, options, sample_info)

        # # Write Execution Record to output directory
        self._write_execution_record(options)

        # If specified, execute the Makefile with nohup
        if options[constants.RUN_SEPARATION]:
            self._run_species_separation(logger, options)

        print(options)
        return(0)




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

    # Write Makefile to output directory
    def _write_makefile(self,logger, options, sample_info):
        print('Wring Rnaseq make file...')
        """
        Write Makefile to results directory to perform species separation.
        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        with fw.writing_to_file(
                fw.MakefileWriter, options[constants.OUTPUT_DIR], "Makefile") as writer:
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
    def _write_execution_record(self,options):
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
            for it in constants.EXECUTION_RECORD_ENTRIES])

        out_file = os.path.join(options[constants.OUTPUT_DIR], "execution_record.txt")
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
        writer.set_variable(constants.NUM_THREADS_VARIABLE, options[constants.NUM_THREADS])
        writer.add_blank_line()

        writer.set_variable(constants.STAR_EXECUTABLE_VARIABLE, options[constants.STAR_EXECUTABLE])
        writer.add_blank_line()

        writer.set_variable(constants.SAMBAMBA_SORT_TMP_DIR_VARIABLE,options[constants.SAMBAMBA_SORT_TMP_DIR])
        writer.add_blank_line()

        sample_names = sample_info.get_sample_names()
        writer.set_variable(constants.SAMPLES_VARIABLE, " ".join(sample_names))
        writer.set_variable(
            constants.RAW_READS_DIRECTORY_VARIABLE,
            options[constants.READS_BASE_DIR] if options[constants.READS_BASE_DIR] else "/")
        writer.set_variable(
            constants.RAW_READS_LEFT_VARIABLE,
            " ".join([",".join(sample_info.get_left_reads(name))
                      for name in sample_names]))

        if sample_info.paired_end_reads():
            writer.set_variable(
                constants.RAW_READS_RIGHT_VARIABLE,
                " ".join([",".join(sample_info.get_right_reads(name))
                          for name in sample_names]))

        writer.add_blank_line()

        for i, species in enumerate(options[constants.SPECIES]):
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
        writer.set_variable(constants.STAR_INDICES_TARGET, "star_indices")
        writer.set_variable(constants.COLLATE_RAW_READS_TARGET, "raw_reads")
        writer.set_variable(constants.MAPPED_READS_TARGET, "mapped_reads")
        writer.set_variable(constants.SORTED_READS_TARGET, "sorted_reads")
        writer.set_variable(constants.FILTERED_READS_TARGET, "filtered_reads")
        writer.add_blank_line()

    def _write_phony_targets(self, logger, writer):
        """
        Write phony target definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                ".PHONY", [constants.ALL_TARGET, constants.CLEAN_TARGET],
                raw_target=True, raw_dependencies=True):
            pass

    def _write_all_target(self, logger, writer):
        """
        Write main target definition to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                constants.ALL_TARGET, [constants.FILTERED_READS_TARGET], raw_target=True):
            pass

    def _write_filtered_reads_target(self, logger, writer, options):
        """
        Write target to separate reads by species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(
                constants.FILTERED_READS_TARGET, [constants.SORTED_READS_TARGET]):
            writer.add_comment(
                "For each sample, take the reads mapping to each genome and " +
                "filter them to their correct species of origin")
            writer.make_target_directory(constants.FILTERED_READS_TARGET)

            writer.add_command("filter_reads", [
                "\"{var}\"".format(var=writer.variable_val(constants.SAMPLES_VARIABLE)),
                writer.variable_val(constants.SORTED_READS_TARGET),
                writer.variable_val(constants.FILTERED_READS_TARGET),
                writer.variable_val(constants.NUM_THREADS_VARIABLE),
                options[constants.MISMATCH_THRESHOLD],
                options[constants.MINMATCH_THRESHOLD],
                options[constants.MULTIMAP_THRESHOLD],
                "--reject-multimaps" if options[constants.REJECT_MULTIMAPS] else "\"\"",
                "{sl}".format(sl=" ".join(options[constants.SPECIES]))])

            if options[constants.DELETE_INTERMEDIATE]:
                writer.remove_target_directory(constants.SORTED_READS_TARGET)


    def _write_sorted_reads_target(self, logger, writer, options):
        """
        Write target to sort reads by name to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        with writer.target_definition(
                constants.SORTED_READS_TARGET, [constants.MAPPED_READS_TARGET]):
            writer.add_comment(
                "For each sample, sort the mapped reads into read name order")
            writer.make_target_directory(constants.SORTED_READS_TARGET)

            writer.add_command(
                "sort_reads",
                ["\"{sl}\"".format(sl=" ".join(options[constants.SPECIES])),
                 "\"{var}\"".format(
                     var=writer.variable_val(constants.SAMPLES_VARIABLE)),
                 writer.variable_val(constants.NUM_THREADS_VARIABLE),
                 writer.variable_val(constants.MAPPED_READS_TARGET),
                 writer.variable_val(constants.SORTED_READS_TARGET),
                 writer.variable_val(constants.SAMBAMBA_SORT_TMP_DIR_VARIABLE)])

            if options[constants.DELETE_INTERMEDIATE]:
                writer.remove_target_directory(constants.MAPPED_READS_TARGET)

    def _write_mapped_reads_target(self, logger, writer, sample_info, options):
        """
        Write target to map reads to each species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """

        index_targets = ["{index}/{species}".format(
            index=writer.variable_val(constants.STAR_INDICES_TARGET),
            species=species)
            for species in options[constants.SPECIES]]

        index_targets.append(writer.variable_val(constants.COLLATE_RAW_READS_TARGET))

        with writer.target_definition(
                constants.MAPPED_READS_TARGET, index_targets, raw_dependencies=True):
            writer.add_comment(
                "Map reads for each sample to each species' genome")
            writer.make_target_directory(constants.MAPPED_READS_TARGET)

            map_reads_params = [
                "\"{sl}\"".format(sl=" ".join(options[constants.SPECIES])),
                "\"{var}\"".format(var=writer.variable_val(constants.SAMPLES_VARIABLE)),
                writer.variable_val(constants.STAR_INDICES_TARGET),
                writer.variable_val(constants.NUM_THREADS_VARIABLE),
                writer.variable_val(constants.COLLATE_RAW_READS_TARGET),
                writer.variable_val(constants.MAPPED_READS_TARGET)]

            map_reads_params.append(
                constants.PAIRED_END_READS_TYPE if sample_info.paired_end_reads()
                else constants.SINGLE_END_READS_TYPE)

            map_reads_params.append(options[constants.STAR_EXECUTABLE])

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
        with writer.target_definition(constants.COLLATE_RAW_READS_TARGET, []):
            writer.add_comment(
                "Create a directory with sub-directories for each sample, " +
                "each of which contains links to the input raw reads files " +
                "for that sample")
            writer.make_target_directory(constants.COLLATE_RAW_READS_TARGET)

            collate_raw_reads_params = [
                "\"{var}\"".format(var=writer.variable_val(constants.SAMPLES_VARIABLE)),
                writer.variable_val(constants.RAW_READS_DIRECTORY_VARIABLE),
                writer.variable_val(constants.COLLATE_RAW_READS_TARGET),
            ]

            if sample_info.paired_end_reads():
                collate_raw_reads_params += [
                    constants.PAIRED_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=writer.variable_val(constants.RAW_READS_LEFT_VARIABLE)),
                    "\"{var}\"".format(
                        var=writer.variable_val(constants.RAW_READS_RIGHT_VARIABLE))
                ]
            else:
                collate_raw_reads_params += [
                    constants.SINGLE_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=writer.variable_val(constants.RAW_READS_LEFT_VARIABLE)),
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
            index=writer.variable_val(constants.STAR_INDICES_TARGET), spec=species)

        with writer.target_definition(target, [], raw_target=True):
            if species_options[constants.GTF_FILE] is None:
                writer.make_target_directory(constants.STAR_INDICES_TARGET)
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
                     writer.variable_val(constants.NUM_THREADS_VARIABLE),
                     target,star_executable])

    def _write_main_star_index_targets(self, logger, writer, options):
        """
        Write targets to create or link to STAR indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        for i, species in enumerate(options[constants.SPECIES]):
            species_options = self._get_species_options(options, i)
            self._write_species_main_star_index_target(
                logger, writer, species, species_options,options[constants.STAR_EXECUTABLE])

    def _write_clean_target(self, logger, writer):
        """
        Write target to clean results directory to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with writer.target_definition(constants.CLEAN_TARGET, [], raw_target=True):
            writer.remove_target_directory(constants.STAR_INDICES_TARGET)
            writer.remove_target_directory(constants.COLLATE_RAW_READS_TARGET)
            writer.remove_target_directory(constants.MAPPED_READS_TARGET)
            writer.remove_target_directory(constants.SORTED_READS_TARGET)
            writer.remove_target_directory(constants.FILTERED_READS_TARGET)


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


    def _write_species_variable_definitions(self, logger, writer, species, species_options):
        """
        Write variable definitions for a particular species.

        logger: logging object
        writer: Makefile writer object
        species: species number string
        species_options: dictionary of options specific to a particular species.
        """
        if species_options[constants.GTF_FILE] is None:
            writer.set_variable(
                self._get_star_index_variable(species), species_options[constants.STAR_INDEX])
        else:
            writer.set_variable(
                self._get_gtf_file_variable(species), species_options[constants.GTF_FILE])
            writer.set_variable(
                self._get_genome_fasta_variable(species), species_options[constants.GENOME_FASTA])

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


class ChipseqSeparator(Separator):
    DOC="""not implemented"""

    # Write Makefile to output directory
    def _write_makefile(self,logger, options, sample_info):
        print('Wring Chipseq make file...')
        # raise NotImplementedError()

    # Write Execution Record to output directory
    def _write_execution_record(self,options):
        print('Wring Rnaseq execution_record...')
        # raise NotImplementedError()


class SeparatorManager(object):
    SEPARATORS = {"rnaseq": RnaseqSeparator,
                  "chipseq": ChipseqSeparator}

    def __init__(self,data_type,args):
        self.args = args
        self.data_type = data_type
        self.separator = self.creat_separator(data_type, args)

    def creat_separator(self,data_type,args):
        separator = self.SEPARATORS[data_type]
        return separator(data_type,args)

    def get_separator(self):
        return self.separator










