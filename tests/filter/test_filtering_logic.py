import sys

from pkg_resources import resource_filename

from sargasso.filter.sample_filter import *


def test_filter_sample_reads_debug(tmpdir):
    read_id = 'SRR7968960.10001094'
    read_id.replace(".", "_")
    sample_name = 'human_1'

    args = ["chipseq",
            "1.0",
            "2.0",
            "999999",
            "mouse",
            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id.replace(".", "_")),
                              "{}___mouse___BLOCK___1.bam".format(sample_name)),
            tmpdir.mkdir("debug").join("{}_mouse_0_filtered.bam".format(sample_name)).strpath,
            "human",
            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id.replace(".", "_")),
                              "{}___human___BLOCK___1.bam".format(sample_name)),
            tmpdir.join("{}_human_0_filtered.bam".format(sample_name)).strpath
            ]

    SampleFilterManager.get(data_type='chipseq').run(args)
    assert 0
    # mouse_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    # rat_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    # human_out = tmpdir.join("chipseq_test_se/chiseq_mouse_se_sample_human_0_filtered.bam").strpath
    # summary_file = tmpdir.join("chipseq_test_se/filtering_result_summary.txt").strpath
    #
    # assert os.path.isfile(mouse_out)
    # assert os.path.isfile(rat_out)
    # assert os.path.isfile(human_out)
    # assert os.path.isfile(summary_file)
    #
    # assert os.access(mouse_out, os.R_OK)
    # assert os.access(rat_out, os.R_OK)
    # assert os.access(human_out, os.R_OK)
    # assert os.access(summary_file, os.R_OK)
    #
    # assert os.stat(mouse_out).st_size != 0
    # assert os.stat(rat_out).st_size != 0
    # assert os.stat(human_out).st_size != 0
    # assert os.stat(summary_file).st_size != 0


def test_filter_sample_reads_batch_chipseq_debug(tmpdir):
    tmpdir.mkdir("debug")
    for type in ['se', 'pe']:
        # for type in ['pe']:
        for sample_name in ['human_1', 'human_2', 'mouse_1', 'mouse_2']:
            # for sample_name in ['human_1']:
            for strategy in ['conservative', 'best', 'recall', 'permissive']:
                # for strategy in ['conservative']:
                #     sample_name='sample_name'
                #     type='pe'
                #     strategy='best'

                read_id = "{}_{}_mis_{}".format(type, strategy, sample_name)
                with open('/home/xinhe/tmp/sargassolog/{}.log'.format(read_id), 'w') as sys.stdout:
                    print('Testing {}'.format(read_id))
                    args = ["chipseq",
                            "1.0",
                            "2.0",
                            "999999",
                            "mouse",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___mouse___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_mouse_0_filtered.bam".format(sample_name)).strpath,
                            "human",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___human___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_human_0_filtered.bam".format(sample_name)).strpath
                            ]

                    SampleFilterManager.get(data_type='chipseq').run(args)


def test_filter_sample_reads_batch_chipseq_without_multimap_debug(tmpdir):
    tmpdir.mkdir("debug")
    for type in ['se', 'pe']:
        # for type in ['pe']:
        for sample_name in ['human_1', 'human_2', 'mouse_1', 'mouse_2']:
            # for sample_name in ['human_1']:
            for strategy in ['conservative', 'best', 'recall', 'permissive']:
                # for strategy in ['conservative']:
                #     sample_name='sample_name'
                #     type='pe'
                #     strategy='best'

                read_id = "{}_{}_mis_{}_without_multimap".format(type, strategy, sample_name)
                with open('/home/xinhe/tmp/sargassolog2/{}.log'.format(read_id), 'w') as sys.stdout:
                    print('Testing {}'.format(read_id))
                    args = ["chipseq",
                            "1.0",
                            "2.0",
                            "999999",
                            "mouse",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___mouse___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_mouse_0_filtered.bam".format(sample_name)).strpath,
                            "human",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___human___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_human_0_filtered.bam".format(sample_name)).strpath
                            ]

                    SampleFilterManager.get(data_type='chipseq').run(args)


def test_filter_sample_reads_batch_rnaseq_debug(tmpdir):
    tmpdir.mkdir("debug")
    for type in ['pe']:
        # for type in ['pe']:
        for sample_name in ['HB1', 'HC1', 'CTR1', 'CTR2']:
            # for sample_name in ['human_1']:
            for strategy in ['conservative', 'best', 'recall', 'permissive']:
                # for strategy in ['conservative']:
                #     sample_name='sample_name'
                #     type='pe'
                #     strategy='best'

                read_id = "{}_{}_mis_{}".format(type, strategy, sample_name)
                with open('/home/xinhe/tmp/sargassolog3/{}.log'.format(read_id), 'w') as sys.stdout:
                    print('Testing {}'.format(read_id))
                    args = ["rnaseq",
                            "1.0",
                            "2.0",
                            "999999",
                            "mouse",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___mouse___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_mouse_0_filtered.bam".format(sample_name)).strpath,
                            "human",
                            resource_filename("tests.data.tmp.{}.filtered_reads.Blocks".format(read_id),
                                              "{}___human___BLOCK___1.bam".format(sample_name)),
                            tmpdir.join("{}_human_0_filtered.bam".format(sample_name)).strpath
                            ]

                    SampleFilterManager.get(data_type='rnaseq').run(args)
