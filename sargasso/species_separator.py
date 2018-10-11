#!/usr/bin/env python
import schema
from commandline_parser import CommandlineParser
from parameter_validator import ParameterValidator
from separators import Separator
from separators import SeparatorManager
from options import Options


def separate_species(args):
    """
    Write a Makefile to perform species separation and optionally execute it.

    args: list of command-line arguments
    """
    # https://github.com/docopt/docopt/blob/master/examples/git/git.py
    data_type = CommandlineParser().parse_datatype(args, Separator.DOC, Options.DATA_TYPE)
    try:
        ParameterValidator.validate_datatype(data_type)
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)

    separator = SeparatorManager().get(data_type)
    separator.run(args)

    # todo remove debug code
    print("--------------------------")
    print("make file content:")
    print("--------------------------")

    with open('/home/xinhe/Projects/Sargasso/results/chipseqtest/Makefile', 'r') as fin:
        print fin.read()
    print("--------------------------")

    print('done')
    exit(1)
