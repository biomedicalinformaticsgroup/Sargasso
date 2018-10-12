import logging
import sys

from factory import Manager
from options import Options


class SargassoLogger(object):

    def __init__(self, logger):
        self.logger = logger

    def init(self, options):
        self._init(sys.stderr, options[Options.LOG_LEVEL_OPTION])

    def _init(self, stream, level):
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
        self.logger.setLevel(Options.LEVELS[level])
        self.logger.addHandler(handler)

    def info(self, msg, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)

    def debug(self, msg, *args, **kwargs):
        self.logger.debug(msg, *args, **kwargs)

    def write_execution_record(self):
        raise NotImplementedError



class LoggerManager(Manager):

    def __init__(self):
        pass

    def get(self):
        logger = logging.getLogger(__name__)
        return SargassoLogger(logger)
