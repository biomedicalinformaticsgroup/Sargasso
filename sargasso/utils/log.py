"""
Utility functions for logging messages. Exports:

get_logger: Return a logger with a specified severity threshold.
"""

import logging
import sys

LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
}

LOG_LEVEL = "log-level"
LOG_LEVEL_OPTION = "--" + LOG_LEVEL


def get_logger(stream, level):
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
    logger.setLevel(LEVELS[level])
    logger.addHandler(handler)

    # fhandler = logging.FileHandler("/home/xinhe/tmp/test.log")
    # fhandler.setFormatter(formatter)
    # logger.addHandler(fhandler)
    return logger


def get_logger_for_options(options):
    """
    Return a Logger instance with a command-line specified severity threshold.

    Return a Logger instance with severity threshold specified by the command
    line option log.LOG_LEVEL. Log messages will be written to standard error.
    options: Dictionary mapping from command-line option strings to option
    values.
    """
    return get_logger(sys.stderr, options[LOG_LEVEL_OPTION])
