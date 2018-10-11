import os
from commandline_parser import CommandlineParserManager
from factory import Manager
from file_writer import MakefileWriterManager
from log import LoggerManager
from options import Options
from parameter_validator import ParameterValidatorManager
from utils import Process


class Separator(object):
    data_type = None

    DOC = """
        Usage: species_separator [-h] <data-type> [<args>...]
            
        options:
           -h
           
        The available commands are:
           rnaseq   rnaseq data
           chipseq  chipseq data
        """

    def __init__(self, data_type):
        self.data_type = data_type

    def run(self, args):
        commandline_parser = CommandlineParserManager().get(self.data_type)
        options = commandline_parser.parse_parameters(args, self.DOC)

        # Validate command-line options
        parameter_validator = ParameterValidatorManager().get(self.data_type)
        parameter_validator.validate(options)

        # Set up logger
        logger = LoggerManager(options).get()

        # Create output directory
        # os.mkdir(options[Options.OUTPUT_DIR])

        # Write Makefile to output directory
        makefile_writer = MakefileWriterManager().get(self.data_type)
        makefile_writer.write(logger, options)

        # Write Execution Record to output directory
        # ExecutionRecordWriter().write(options)

        # If specified, execute the Makefile with nohup
        if options[Options.RUN_SEPARATION]:
            self._run_species_separation(options)

    def _run_species_separation(self, options):
        """
        Executes the written Makefile with nohup.

        logger: logging object
        options: dictionary of command-line options
        """
        Process.run_in_directory(options[Options.OUTPUT_DIR], "make")


class RnaseqSeparator(Separator):
    data_type = 'rnaseq'
    DOC = """
Usage:
    species_separator <data-type> 
    [--log-level=<log-level>]
    [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
    [--mismatch-threshold=<mismatch-threshold>]
    [--minmatch-threshold=<minmatch-threshold>]
    [--multimap-threshold=<multimap-threshold>]
    [--reject-multimaps]
    [--best] [--conservative] [--recall] [--permissive]
    [--run-separation]
    [--delete-intermediate]
    [--mapper-executable=<mapper-executable>]
    [--mapper-index-executable=<mapper-executable>]
    [--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>]
    <samples-file> <output-dir>
    (<species> <species-info>)
    (<species> <species-info>)
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
    <species-info>
        Either a STAR index directory for the species, or a comma-separated list of
        (i) a GTF annotation file and (ii) a directory containing genome FASTA
        files for the species.
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
    --mapper-executable=<mapper-executable>
        Specify STAR executable path. Use this to run Sargasso with a particular
        version of STAR [default: STAR].
    --mapper-index-executable=<mapper-index-executable>
        same as <mapper-executable>  [default: STAR].
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
    
    species_separator --reads-base-dir=/srv/data/rnaseq --num-threads 4 \
                        --run-separation samples.tsv my_results \
                        mouse /srv/data/genome/mouse/STAR_Index \
                        rat /srv/data/genome/rat/STAR_Index
    
"""


class ChipseqSeparator(Separator):
    data_type = 'chipseq'
    DOC = """
Usage:
    species_separator <data-type> 
    [--log-level=<log-level>]
    [--reads-base-dir=<reads-base-dir>] [--num-threads=<num-threads>]
    [--mismatch-threshold=<mismatch-threshold>]
    [--minmatch-threshold=<minmatch-threshold>]
    [--multimap-threshold=<multimap-threshold>]
    [--reject-multimaps]
    [--best] [--conservative] [--recall] [--permissive]
    [--run-separation]
    [--delete-intermediate]
    [--mapper-executable=<mapper-executable>]
    [--mapper-index-executable=<mapper-index-executable>]
    [--sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>]
    <samples-file> <output-dir>
    (<species> <species-info>)
    (<species> <species-info>)
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
    <species-info>
        Either a bowtie2 index directory for the species, or a the genome FASTA
        files for the species.
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
    --mapper-executable=<mapper-executable>
        Specify bowtie2 executable path. Use this to run Sargasso with a particular
        version of bowtie2 [default: bowtie2].
    --mapper-index-executable=<mapper-executable>
        Specify bowtie2 executable path. Use this to run Sargasso with a particular
        version of bowtie2 [default: bowtie2-build].
    --sambamba-sort-tmp-dir=<sambamba-sort-tmp-dir>
        Specify 'sambamba sort' temporary folder path [default: /tmp].
    
    species_separator chipseq --reads-base-dir=/home/xinhe/Projects/Sargasso/pipeline_test/data/fastq/ --num-threads 4 \
                        --run-separation /home/xinhe/Projects/Sargasso/pipeline_test/data/fastq/chipseq.tsv \
                        my_results \
                        mouse /srv/data/genome/mouse/ensembl-93/bowtie2_indexes/primary_assembly \
                        human /srv/data/genome/human/ensembl-93/bowtie2_indexes/primary_assembly \
                        rat /srv/data/genome/rat/ensembl-93/bowtie2_indexes/toplevel
    
"""


class SeparatorManager(Manager):
    SEPARATORS = {"rnaseq": RnaseqSeparator,
                  "chipseq": ChipseqSeparator}

    def get(self, data_type):
        return self.SEPARATORS[data_type](data_type)
