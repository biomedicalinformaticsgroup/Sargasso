import docopt
import textwrap




from .__init__ import __version__
from options import Options
from factory import Manager


class CommandlineParser(object):
    def __init__(self):
        # self.DATA_TYPE=data_type
        # self.args = args
        # self.doc = doc
        # self.pv = ParameterValidator()
        pass

    @staticmethod
    def parse(args, doc, options_first=False):
        # Read in command-line options
        docstring = CommandlineParser.substitute_common_options_into_usage(doc)
        options = docopt.docopt(docstring, argv=args, version="species_separator v" + __version__,options_first=options_first)
        return options


    @staticmethod
    def substitute_common_options_into_usage(usage_msg, **substitutions):
        """
        Substitute common option and other interpolations into a usage message.

        Substitute help, version and logging level option specifications and
        descriptions into a script's usage message; also substitute any other
        interpolations specified by additional keyword arguments.
        usage_msg: A script's usage message.
        substitutions: Additional key=value interpolations to substitute into the
        usage message.
        """
        help_spec = "-h --help"
        help_desc = "Show this message."

        ver_spec = "-v --version"
        ver_desc = "Show version."

        log_spec = "{log_option}=<{log_level}>".format(
            log_option=Options.LOG_LEVEL_OPTION, log_level=Options.LOG_LEVEL)
        log_desc = ("Set logging level " +
                    "(one of {log_level_vals}) [default: info].").format(
            log_level_vals=str(Options.LEVELS.keys()))
        log_desc = '\n'.join(textwrap.wrap(
            log_desc, width=75, subsequent_indent="    "))

        return usage_msg.format(
            log_option_spec=log_spec, log_option_description=log_desc,
            help_option_spec=help_spec, help_option_description=help_desc,
            ver_option_spec=ver_spec, ver_option_description=ver_desc,
            **substitutions)


    def get_sample_info(self):
        raise NotImplementedError
        return sample_info

    def get_logger(self):
        raise NotImplementedError
        return self.logger
class RnaseqCommandlineParser(CommandlineParser):
    pass
class ChipseCommandlineParser(CommandlineParser):
    pass

class CommandlineParserManager(Manager):
    PARSERS = {"rnaseq" : RnaseqCommandlineParser,
               "chipseq": ChipseCommandlineParser}

    def __init__(self):
        # self.args = args
        # self.data_type = data_type
        # self.doc = doc
        # self.parser=self._create()
        pass

    def _create(self,data_type):
        commandlineParser = self.PARSERS[data_type]
        return commandlineParser()

    def get(self,data_type):
        return self._create(data_type)
