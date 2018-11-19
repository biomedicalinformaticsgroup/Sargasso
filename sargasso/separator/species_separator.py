#!/usr/bin/env python
import schema

from sargasso.separator.commandline_parser import CommandlineParser
from sargasso.separator.options import Options
from sargasso.separator.parameter_validator import ParameterValidator
from sargasso.separator.separators import Separator
from sargasso.separator.separators import SeparatorManager


def separate_species(args):
    """
    Write a Makefile to perform species separation and optionally execute it.

    args: list of command-line arguments
    """
    # https://github.com/docopt/docopt/blob/master/examples/git/git.py
    data_type = CommandlineParser.parse_datatype(args, Separator.DOC, Options.DATA_TYPE)
    try:
        ParameterValidator.validate_datatype(data_type)
    except schema.SchemaError as exc:
        exit("Exiting: " + exc.code)

    separator = SeparatorManager.get(data_type)
    separator.run(args)
