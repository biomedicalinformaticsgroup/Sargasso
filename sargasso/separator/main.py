from sargasso.filter.filter_controllers import FilterController
from sargasso.filter.sample_filterer import SampleFilter
from sargasso.separator.data_types import get_data_type_manager
from sargasso.separator.separators import Separator


def separate_species(args):
    data_type_manager = get_data_type_manager(args, Separator.DOC)
    data_type_manager.get_separator().run(args)


def filter_control(args):
    data_type_manager = get_data_type_manager(args, FilterController.DOC)
    data_type_manager.get_filter_controller().run(args)


def filter_sample_reads(args):
    data_type_manager = get_data_type_manager(args, SampleFilter.DOC)
    data_type_manager.get_sample_filterer().run(args)
