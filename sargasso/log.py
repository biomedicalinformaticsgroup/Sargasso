import sys
import logging
from factory import Manager
from options import Options


class Logger(object):

    def write_execution_record(self):
        raise NotImplementedError


class LoggerManager(Manager):

    def __init__(self, options):
        self.options = options
        self.logger = self._create(options)

    def _create(self, options):
        """
        Return a Logger instance with a command-line specified severity threshold.

        Return a Logger instance with severity threshold specified by the command
        line option log.LOG_LEVEL. Log messages will be written to standard out.
        options: Dictionary mapping from command-line option strings to option
        values.
        """
        return self._get_logger(sys.stderr, options[Options.LOG_LEVEL_OPTION])

    def _get_logger(self, stream, level):
        """
        Return a Logger instance with the specified severity threshold.

        Return a Logger instance with the specified severity threshold, where the
        threshold level should be a key of the 'LEVELS' dictionary. Log messages
        will contain the current time and message severity level.
        stream: Output stream to which the logger will write messages.
        level: Severity threshold level, which should be a key of the 'LEVELS'
        dictionary.
        """
        formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s: %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')
        handler = logging.StreamHandler(stream)
        handler.setFormatter(formatter)
        logger = logging.getLogger(__name__)
        logger.setLevel(Options.LEVELS[level])
        logger.addHandler(handler)
        return logger

    def get(self):
        return self.logger
