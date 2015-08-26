#!/usr/bin/env python

"""Usage:
    filter_mapped_hits.py [--log-level=<log-level>] <species-one> <species-one-input-bam> <species-one-output-bam> <species-two> <species-two-input-bam> <species-two-output-bam>

-h --help                   Show this message.
-v --version                Show version.
--log-level=<log-level>     Set logging level (one of {log_level_vals}) [default: info].
<species-one>               Name of first species.
<species-one-input-bam>     BAM file containing read hits against first species' genome.
<species-one-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
<species-two>               Name of first species.
<species-two-input-bam>     BAM file containing read hits against first species' genome.
<species-two-output-bam>    BAM file to which reads assigned to first species after filtering will be written.
"""

import docopt
import log

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LEVELS.keys())


def filter_mapped_hits(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring, version="filter_mapped_hits v0.1")
    print(options)
