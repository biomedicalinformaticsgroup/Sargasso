import os
import os.path
import schema
from . import parameter_defination as pd


from . import file_writer as fw
from datetime import datetime

class makefile_writer(object):


    def __init__(self):
        pass
        # self.logger = logger
        # self.options = options
        # self.sample_info = sample_info

    def write_makefile(self, logger, options, sample_info):
        raise NotImplementedError()






class rnaseq_makefile_writer(makefile_writer):

    def write_makefile(self, logger, options, sample_info):
        """
        Write Makefile to results directory to perform species separation.
        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        with fw.writing_to_file(
                fw.MakefileWriter, options[pd.OUTPUT_DIR], "Makefile") as writer:
            self._write_variable_definitions(logger, writer, options, sample_info)
            self._write_target_variable_definitions(logger, writer)
            self._write_phony_targets(logger, writer)
            self._write_all_target(logger, writer)
            self._write_filtered_reads_target(logger, writer, options)
            self._write_sorted_reads_target(logger, writer, options)
            self._write_mapped_reads_target(logger, writer, sample_info, options)
            # self._write_masked_reads_target(logger, writer)
            self._write_collate_raw_reads_target(logger, writer, sample_info)
            # self._write_mask_star_index_targets(logger, writer, options)
            self._write_main_star_index_targets(logger, writer, options)
            self._write_clean_target(logger, writer)

    def _get_species_options(self, species_index):
        """
        Return a dictionary containing command-line options for a species.

        options: dictionary containing all command-line options.
        species_index: which species to return options for
        """

        species_options = { self.SPECIES_NAME: self.options[self.SPECIES][species_index] }

        star_infos = self.options[self.SPECIES_STAR_INFO][species_index].split(",")

        if len(star_infos) == 1:
            species_options[self.STAR_INDEX] = star_infos[0]
            species_options[self.GTF_FILE] = None
            species_options[self.GENOME_FASTA] = None
        elif len(star_infos) == 2:
            species_options[self.STAR_INDEX] = None
            species_options[self.GTF_FILE] = star_infos[0]
            species_options[self.GENOME_FASTA] = star_infos[1]
        else:
            raise schema.SchemaError(
                None, "Should specify either a STAR index or both GTF file " +
                      "and genome FASTA directory for species {species}.".
                      format(species=species_options[self.SPECIES_NAME]))

        return species_options

    def _get_gtf_file_variable(species):
        """
        Return GTF file variable name for a species.

        species: species identifier
        """
        return "{species}_GTF".format(species=species.upper())


    def _get_genome_fasta_variable(species):
        """
        Return genome FASTA variable name for a species.

        species: species identifier
        """
        return "{species}_GENOME_FASTA_DIR".format(species=species.upper())


    def _get_star_index_variable(species):
        """
        Return STAR index variable name for a species.

        species: species identifier
        """
        return "{species}_STAR_DIR".format(species=species.upper())


    def _write_variable_definitions(self):
        """
        Write variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        self.writer.set_variable(self.NUM_THREADS_VARIABLE, self.options[self.NUM_THREADS])
        self.writer.add_blank_line()

        self.writer.set_variable(self.STAR_EXECUTABLE_VARIABLE, self.options[self.STAR_EXECUTABLE])
        self.writer.add_blank_line()

        self.writer.set_variable(self.SAMBAMBA_SORT_TMP_DIR_VARIABLE,self.options[self.SAMBAMBA_SORT_TMP_DIR])
        self.writer.add_blank_line()

        sample_names = self.sample_info.get_sample_names()
        self.writer.set_variable(self.SAMPLES_VARIABLE, " ".join(sample_names))
        self.writer.set_variable(
            self.RAW_READS_DIRECTORY_VARIABLE,
            self.options[self.READS_BASE_DIR] if self.options[self.READS_BASE_DIR] else "/")
        self.writer.set_variable(
            self.RAW_READS_LEFT_VARIABLE,
            " ".join([",".join(self.sample_info.get_left_reads(name))
                      for name in sample_names]))

        if self.sample_info.paired_end_reads():
            self.writer.set_variable(
                self.RAW_READS_RIGHT_VARIABLE,
                " ".join([",".join(self.sample_info.get_right_reads(name))
                          for name in sample_names]))

        self.writer.add_blank_line()

        for i, species in enumerate(self.options[self.SPECIES]):
            species_options = self._get_species_options(i)
            self._write_species_variable_definitions(
                self.logger, self.writer, species, species_options)

        self.writer.add_blank_line()


    def _write_target_variable_definitions(self):
        """
        Write target directory variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        self.writer.set_variable(self.STAR_INDICES_TARGET, "star_indices")
        self.writer.set_variable(self.COLLATE_RAW_READS_TARGET, "raw_reads")
        self.writer.set_variable(self.MAPPED_READS_TARGET, "mapped_reads")
        self.writer.set_variable(self.SORTED_READS_TARGET, "sorted_reads")
        self.writer.set_variable(self.FILTERED_READS_TARGET, "filtered_reads")
        self.writer.add_blank_line()


    def _write_phony_targets(self):
        """
        Write phony target definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.writer.target_definition(
                ".PHONY", [self.ALL_TARGET, self.CLEAN_TARGET],
                raw_target=True, raw_dependencies=True):
            pass


    def _write_all_target(self):
        """
        Write main target definition to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.writer.target_definition(
                self.ALL_TARGET, [self.FILTERED_READS_TARGET], raw_target=True):
            pass


    def _write_filtered_reads_target(self):
        """
        Write target to separate reads by species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.writer.target_definition(
                self.FILTERED_READS_TARGET, [self.SORTED_READS_TARGET]):
            self.writer.add_comment(
                "For each sample, take the reads mapping to each genome and " +
                "filter them to their correct species of origin")
            self.writer.make_target_directory(self.FILTERED_READS_TARGET)

            self.writer.add_command("filter_reads", [
                "\"{var}\"".format(var=self.writer.variable_val(self.SAMPLES_VARIABLE)),
                self.writer.variable_val(self.SORTED_READS_TARGET),
                self.writer.variable_val(self.FILTERED_READS_TARGET),
                self.writer.variable_val(self.NUM_THREADS_VARIABLE),
                self.options[self.MISMATCH_THRESHOLD],
                self.options[self.MINMATCH_THRESHOLD],
                self.options[self.MULTIMAP_THRESHOLD],
                "--reject-multimaps" if self.options[self.REJECT_MULTIMAPS] else "\"\"",
                "{sl}".format(sl=" ".join(self.options[self.SPECIES]))])

            if self.options[self.DELETE_INTERMEDIATE]:
                self.writer.remove_target_directory(self.SORTED_READS_TARGET)


    def _write_sorted_reads_target(self):
        """
        Write target to sort reads by name to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        with self.writer.target_definition(
                self.SORTED_READS_TARGET, [self.MAPPED_READS_TARGET]):
            self.writer.add_comment(
                "For each sample, sort the mapped reads into read name order")
            self.writer.make_target_directory(self.SORTED_READS_TARGET)

            self.writer.add_command(
                "sort_reads",
                ["\"{sl}\"".format(sl=" ".join(self.options[self.SPECIES])),
                 "\"{var}\"".format(
                     var=self.writer.variable_val(self.SAMPLES_VARIABLE)),
                 self.writer.variable_val(self.NUM_THREADS_VARIABLE),
                 self.writer.variable_val(self.MAPPED_READS_TARGET),
                 self.writer.variable_val(self.SORTED_READS_TARGET),
                 self.writer.variable_val(self.SAMBAMBA_SORT_TMP_DIR_VARIABLE)])

            if self.options[self.DELETE_INTERMEDIATE]:
                self.writer.remove_target_directory(self.MAPPED_READS_TARGET)


    def _write_mapped_reads_target(self):
        """
        Write target to map reads to each species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """

        index_targets = ["{index}/{species}".format(
            index=self.writer.variable_val(self.STAR_INDICES_TARGET),
            species=species)
            for species in self.options[self.SPECIES]]

        index_targets.append(self.writer.variable_val(self.COLLATE_RAW_READS_TARGET))

        with self.writer.target_definition(
                self.MAPPED_READS_TARGET, index_targets, raw_dependencies=True):
            self.writer.add_comment(
                "Map reads for each sample to each species' genome")
            self.writer.make_target_directory(self.MAPPED_READS_TARGET)

            map_reads_params = [
                "\"{sl}\"".format(sl=" ".join(self.options[self.SPECIES])),
                "\"{var}\"".format(var=self.writer.variable_val(self.SAMPLES_VARIABLE)),
                self.writer.variable_val(self.STAR_INDICES_TARGET),
                self.writer.variable_val(self.NUM_THREADS_VARIABLE),
                self.writer.variable_val(self.COLLATE_RAW_READS_TARGET),
                self.writer.variable_val(self.MAPPED_READS_TARGET)]

            map_reads_params.append(
                self.PAIRED_END_READS_TYPE if self.sample_info.paired_end_reads()
                else self.SINGLE_END_READS_TYPE)

            map_reads_params.append(self.options[self.STAR_EXECUTABLE])

            self.writer.add_command("map_reads", map_reads_params)


            #def _write_masked_reads_target(logger, writer, options):
            # TODO: Write "masked_reads" target
            #pass


    def _write_collate_raw_reads_target(self):
        """
        Write target to collect raw reads files to Makefile.

        logger: logging object
        writer: Makefile writer object
        sample_info: object encapsulating samples and their accompanying read files
        """
        with self.writer.target_definition(self.COLLATE_RAW_READS_TARGET, []):
            self.writer.add_comment(
                "Create a directory with sub-directories for each sample, " +
                "each of which contains links to the input raw reads files " +
                "for that sample")
            self.writer.make_target_directory(self.COLLATE_RAW_READS_TARGET)

            collate_raw_reads_params = [
                "\"{var}\"".format(var=self.writer.variable_val(self.SAMPLES_VARIABLE)),
                self.writer.variable_val(self.RAW_READS_DIRECTORY_VARIABLE),
                self.writer.variable_val(self.COLLATE_RAW_READS_TARGET),
            ]

            if self.sample_info.paired_end_reads():
                collate_raw_reads_params += [
                    self.PAIRED_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=self.writer.variable_val(self.RAW_READS_LEFT_VARIABLE)),
                    "\"{var}\"".format(
                        var=self.writer.variable_val(self.RAW_READS_RIGHT_VARIABLE))
                ]
            else:
                collate_raw_reads_params += [
                    self.SINGLE_END_READS_TYPE,
                    "\"{var}\"".format(
                        var=self.writer.variable_val(self.RAW_READS_LEFT_VARIABLE)),
                    "\"\""
                ]

            self.writer.add_command("collate_raw_reads", collate_raw_reads_params)


            #def _write_mask_star_index_targets(logger, writer, options):
            # TODO: Write targets for STAR indices for mask sequences
            #pass


    def _write_species_main_star_index_target(self, species, species_options):
        """
        Write target to create or link to STAR index for a species to Makefile.

        logger: logging object
        writer: Makefile writer object
        species_var: Makefile variable for species
        species_options: dictionary of command-line options for the particular
        species
        """
        target = "{index}/{spec}".format(
            index=self.writer.variable_val(self.STAR_INDICES_TARGET), spec=species)

        with self.writer.target_definition(target, [], raw_target=True):
            if species_options[self.GTF_FILE] is None:
                self.writer.make_target_directory(self.STAR_INDICES_TARGET)
                self.writer.add_command(
                    "ln",
                    ["-s",
                     self.writer.variable_val(self._get_star_index_variable(species)),
                     target])
            else:
                self.writer.make_target_directory(target, raw_target=True)
                self.writer.add_command(
                    "build_star_index",
                    [self.writer.variable_val(self._get_genome_fasta_variable(species)),
                     self.writer.variable_val(self._get_gtf_file_variable(species)),
                     self.writer.variable_val(self.NUM_THREADS_VARIABLE),
                     target,self.options[self.STAR_EXECUTABLE]])


    def _write_main_star_index_targets(self):
        """
        Write targets to create or link to STAR indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        for i, species in enumerate(self.options[self.SPECIES]):
            species_options = self._get_species_options(i)
            self._write_species_main_star_index_target( species, species_options)


    def _write_clean_target(self):
        """
        Write target to clean results directory to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        with self.writer.target_definition(self.CLEAN_TARGET, [], raw_target=True):
            self.writer.remove_target_directory(self.STAR_INDICES_TARGET)
            self.writer.remove_target_directory(self.COLLATE_RAW_READS_TARGET)
            self.writer.remove_target_directory(self.MAPPED_READS_TARGET)
            self.writer.remove_target_directory(self.SORTED_READS_TARGET)
            self.writer.remove_target_directory(self.FILTERED_READS_TARGET)


    def write_makefile(self):
        """
        Write Makefile to results directory to perform species separation.

        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        with fw.writing_to_file(
                fw.MakefileWriter, self.options[self.OUTPUT_DIR], "Makefile") as writer:
            self._write_variable_definitions()
            self._write_target_variable_definitions()
            self._write_phony_targets()
            self._write_all_target()
            self._write_filtered_reads_target()
            self._write_sorted_reads_target()
            self._write_mapped_reads_target()
            #_write_masked_reads_target(logger, writer)
            self._write_collate_raw_reads_target()
            #_write_mask_star_index_targets(logger, writer, options)
            self._write_main_star_index_targets()
            self._write_clean_target()


    def write_execution_record(self):
        """
        Write a log file containing all execution parameters in addition to the
        time and date of execution

        options: dictionary of command-line options
        """
        out_text = "Execution Record - {t}\n".format(
            t=str(datetime.now().isoformat()))
        out_text += "\n".join(["{desc}: {val}".format(
            desc=it[0], val=str(self.options[it[1]]))
            for it in self.EXECUTION_RECORD_ENTRIES])

        out_file = os.path.join(self.options[self.OUTPUT_DIR], "execution_record.txt")
        with open(out_file, 'w') as erf:
            erf.write(out_text)


class chipseq_makefile_writer(makefile_writer):


    def _write_variable_definitions(self):
        """
        Write variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        raise NotImplementedError()


    def _write_target_variable_definitions(logger, writer):
        """
        Write target directory variable definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        raise NotImplementedError()


    def _write_phony_targets(logger, writer):
        """
        Write phony target definitions to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        raise NotImplementedError()

    def _write_all_target(logger, writer):
        """
        Write main target definition to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        raise NotImplementedError()


    def _write_filtered_reads_target(logger, writer, options):
        """
        Write target to separate reads by species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        raise NotImplementedError()


    def _write_sorted_reads_target(logger, writer, options):
        """
        Write target to sort reads by name to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        raise NotImplementedError()

    def _write_mapped_reads_target(logger, writer, sample_info, options):
        """
        Write target to map reads to each species to Makefile.

        logger: logging object
        writer: Makefile writer object
        """

        raise NotImplementedError()


    def _write_collate_raw_reads_target(logger, writer, sample_info):
        """
        Write target to collect raw reads files to Makefile.

        logger: logging object
        writer: Makefile writer object
        sample_info: object encapsulating samples and their accompanying read files
        """
        raise NotImplementedError()


    def _write_species_main_star_index_target(
            logger, writer, species, species_options,star_executable):
        """
        Write target to create or link to STAR index for a species to Makefile.

        logger: logging object
        writer: Makefile writer object
        species_var: Makefile variable for species
        species_options: dictionary of command-line options for the particular
        species
        """
        raise NotImplementedError()


    def _write_main_star_index_targets(logger, writer, options):
        """
        Write targets to create or link to STAR indices to Makefile.

        logger: logging object
        writer: Makefile writer object
        options: dictionary of command-line options
        """
        raise NotImplementedError()


    def _write_clean_target(logger, writer):
        """
        Write target to clean results directory to Makefile.

        logger: logging object
        writer: Makefile writer object
        """
        raise NotImplementedError()


    def write_makefile(self):
        """
        Write Makefile to results directory to perform species separation.

        logger: logging object
        options: dictionary of command-line options
        sample_info: object encapsulating samples and their accompanying read files
        """
        raise NotImplementedError()


    def write_execution_record(options):
        """
        Write a log file containing all execution parameters in addition to the
        time and date of execution

        options: dictionary of command-line options
        """
        raise NotImplementedError()