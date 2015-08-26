#!/usr/bin/env python

"""Usage:
    filter_mapped_hits.py [--log-level=<log-level>] <species-one> <species-one-input-bam> <species-one-output-bam> <species-two> <species-two-input-bam> <species-two-output-bam>

Options:
{help_option_spec}
    {help_option_description}
{ver_option_spec}
    {ver_option_description}
{log_option_spec}
    {log_option_description}
<species-one>               Name of first species.
<species-one-input-bam>     BAM file containing read hits against first species' genome.
<species-one-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<species-two>               Name of first species.
<species-two-input-bam>     BAM file containing read hits against first species' genome.
<species-two-output-bam>    BAM file to which reads assigned to first species after filtering will be written.

TODO: what does this script do...
"""

import docopt

from . import options as opt
from .__init__ import __version__


def _validate_command_line_options(options):
    # TODO: validate log level
    # TODO: check both input BAM files exist
    pass


def _filter_mapped_hits(logger, options):
    # 1. Attempt to read all hits for the first/next read in each species'
    # input BAM file. Go to (2).
    # 2. If no more reads can be read for species 1:
    #       - all remaining reads in the input file for species 2 can be written
    #       to the output file for species 2, or discarded as inherently
    #       potentially ambiguous.
    #       - END.
    #    Else go to (3).
    # 3. If no more reads can be read for species 2:
    #       - all remaining reads in the input file for species 1 can be
    #       written to the output file for species 1, or discarded as
    #       inherently potentially ambiguous.
    #       - END.
    #   Else go to (4).
    # 4. We have the hits for a read for each species. Examine the names of the
    # reads in species 1 and species 2. Go to (5).
    # 5. If the read names are the same:
    #       - compare the hits for each species to determine which species to
    #       assign the read to.
    #       - write the hits to the output file for the correct species (or
    #       neither, if the read cannot be unambiguously assigned).
    #       - Go to (1).
    #   Else go to (6).
    # 6. If species 1 read name < species 2 read name, there are no hits for
    # the species 1 read in species 2:
    #       - write the hits for this read to the output file for species 1, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 2. If no more
    #       reads can be read for species 2, go to (2). Else go to (4).
    #   Else go to (7)
    # 7. If species 2 read name < species 1 read name, there are no hits for
    # the speices 2 read in species 1:
    #       - write the hits for this read to the output file for species 2, or
    #       discard as inherently potentially ambiguous.
    #       - attempt to read hits for the next read for species 1. If no more
    #       reads can be read for species 1, go to (3). Else go to (4).
    # n.b. throughout, we should always keep a record of the last read name
    # read for a species; when hits for the next read are read for that
    # species, we should check the reads are in name order, and exit with a
    # failure code if not.
    pass


def filter_mapped_hits(args):
    # Read in command-line options
    docstring = opt.substitute_common_options_into_usage(__doc__)
    options = docopt.docopt(docstring, argv=args,
                            version="filter_mapped_hits v" + __version__)

    # Validate command-line options
    _validate_command_line_options(options)

    # Set up logger
    logger = opt.get_logger_for_options(options)

    _filter_mapped_hits(logger, options)
