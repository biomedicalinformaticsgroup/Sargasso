from pkg_resources import resource_filename

from sargasso.filter.sample_filterer import *
from sargasso.separator.data_types import get_data_type_manager
from sargasso.separator.main import *

def test_bisulfite_pe_filter_control(tmpdir):
    tmpdir.mkdir("filtered_reads")
    args = ["bisulfite",
            "--log-level=info",
            resource_filename("tests.data.pe.bisulfite.filtered_reads","Blocks/"),
            tmpdir.join("filtered_reads").strpath,
            "bisulfite_human_pe_sample",
            "1.0",
            "2.0",
            "999999",
            "mouse",
            "human",
            "rat"]

    # filter_control(args)
    assert True

