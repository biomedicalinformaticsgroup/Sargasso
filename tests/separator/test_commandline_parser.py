from sargasso.separator.commandline_parser import *
from sargasso.separator.separators import Separator


def test_parse_datatype():
    data_type = CommandlineParser.parse_datatype(["rnaseq", "options"], Separator.DOC, '<data-type>',
                                                 options_first=True)
    assert data_type == "rnaseq"

    data_type = CommandlineParser.parse_datatype(["chipseq", "options"], Separator.DOC, '<data-type>',
                                                 options_first=True)
    assert data_type == "chipseq"
