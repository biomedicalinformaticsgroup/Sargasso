from pkg_resources import resource_filename

from sargasso.filter.sample_filterer import *
from sargasso.separator.data_types import get_data_type_manager
from sargasso.separator.separators import *


def test_filter_sample_reads_pe_rnaseq(tmpdir):
    args = ["rnaseq",
            "1.0",
            "2.0",
            "999999",
            "--log-level=info",
            "mouse",
            resource_filename("tests.data.pe.rnaseq.filtered_reads.Blocks",
                              "rnaseq_mouse_rat_sample___mouse___BLOCK___1.bam"),
            tmpdir.mkdir("rnaseq_test").join("rnaseq_mouse_rat_sample_mouse_0_filtered.bam").strpath,
            "human",
            resource_filename("tests.data.pe.rnaseq.filtered_reads.Blocks",
                              "rnaseq_mouse_rat_sample___human___BLOCK___1.bam"),
            tmpdir.join("rnaseq_test/rnaseq_mouse_rat_sample_human_0_filtered.bam").strpath,
            "rat",
            resource_filename("tests.data.pe.rnaseq.filtered_reads.Blocks",
                              "rnaseq_mouse_rat_sample___rat___BLOCK___1.bam"),
            tmpdir.join("rnaseq_test/rnaseq_mouse_rat_sample_rat_0_filtered.bam").strpath
            ]

    get_data_type_manager(['rnaseq'], SampleFilterer.DOC).get_sample_filterer().run(args)

    mouse_out = tmpdir.join("rnaseq_test/rnaseq_mouse_rat_sample_human_0_filtered.bam").strpath
    rat_out = tmpdir.join("rnaseq_test/rnaseq_mouse_rat_sample_human_0_filtered.bam").strpath
    human_out = tmpdir.join("rnaseq_test/rnaseq_mouse_rat_sample_human_0_filtered.bam").strpath
    summary_file = tmpdir.join("rnaseq_test/filtering_result_summary.txt").strpath

    assert os.path.isfile(mouse_out)
    assert os.path.isfile(rat_out)
    assert os.path.isfile(human_out)
    assert os.path.isfile(summary_file)

    assert os.access(mouse_out, os.R_OK)
    assert os.access(rat_out, os.R_OK)
    assert os.access(human_out, os.R_OK)
    assert os.access(summary_file, os.R_OK)

    assert os.stat(mouse_out).st_size != 0
    assert os.stat(rat_out).st_size != 0
    assert os.stat(human_out).st_size != 0
    assert os.stat(summary_file).st_size != 0


def test_filter_sample_reads_pe_chipseq(tmpdir):
    args = ["dnaseq",
            "1.0",
            "2.0",
            "999999",
            "mouse",
            resource_filename("tests.data.pe.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_sample___mouse___BLOCK___1.bam"),
            tmpdir.mkdir("chipseq_test_pe").join("chiseq_mouse_sample_mouse_0_filtered.bam").strpath,
            "human",
            resource_filename("tests.data.pe.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_sample___human___BLOCK___1.bam"),
            tmpdir.join("chipseq_test_pe/chiseq_mouse_sample_human_0_filtered.bam").strpath,
            "rat",
            resource_filename("tests.data.pe.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_sample___rat___BLOCK___1.bam"),
            tmpdir.join("chipseq_test_pe/chiseq_mouse_sample_rat_0_filtered.bam").strpath
            ]

    get_data_type_manager(['dnaseq'], SampleFilterer.DOC).get_sample_filterer().run(args)

    mouse_out = tmpdir.join("chipseq_test_pe/chiseq_mouse_sample_human_0_filtered.bam").strpath
    rat_out = tmpdir.join("chipseq_test_pe/chiseq_mouse_sample_human_0_filtered.bam").strpath
    human_out = tmpdir.join("chipseq_test_pe/chiseq_mouse_sample_human_0_filtered.bam").strpath
    summary_file = tmpdir.join("chipseq_test_pe/filtering_result_summary.txt").strpath

    assert os.path.isfile(mouse_out)
    assert os.path.isfile(rat_out)
    assert os.path.isfile(human_out)
    assert os.path.isfile(summary_file)

    assert os.access(mouse_out, os.R_OK)
    assert os.access(rat_out, os.R_OK)
    assert os.access(human_out, os.R_OK)
    assert os.access(summary_file, os.R_OK)

    assert os.stat(mouse_out).st_size != 0
    assert os.stat(rat_out).st_size != 0
    assert os.stat(human_out).st_size != 0
    assert os.stat(summary_file).st_size != 0


def test_filter_sample_reads_se_chipseq(tmpdir):
    args = ["dnaseq",
            "1.0",
            "2.0",
            "999999",
            "mouse",
            resource_filename("tests.data.se.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_se_sample___mouse___BLOCK___1.bam"),
            tmpdir.mkdir("chipseq_test_se").join("chiseq_mouse_se_sample_mouse_0_filtered.bam").strpath,
            "human",
            resource_filename("tests.data.se.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_se_sample___human___BLOCK___1.bam"),
            tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath,
            "rat",
            resource_filename("tests.data.se.chipseq.filtered_reads.Blocks",
                              "chiseq_mouse_se_sample___rat___BLOCK___1.bam"),
            tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_rat_0_filtered.bam").strpath
            ]

    get_data_type_manager(['dnaseq'], SampleFilterer.DOC).get_sample_filterer().run(args)

    mouse_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    rat_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    human_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    summary_file = tmpdir.join("chipseq_test_se/filtering_result_summary.txt").strpath

    assert os.path.isfile(mouse_out)
    assert os.path.isfile(rat_out)
    assert os.path.isfile(human_out)
    assert os.path.isfile(summary_file)

    assert os.access(mouse_out, os.R_OK)
    assert os.access(rat_out, os.R_OK)
    assert os.access(human_out, os.R_OK)
    assert os.access(summary_file, os.R_OK)

    assert os.stat(mouse_out).st_size != 0
    assert os.stat(rat_out).st_size != 0
    assert os.stat(human_out).st_size != 0
    assert os.stat(summary_file).st_size != 0
