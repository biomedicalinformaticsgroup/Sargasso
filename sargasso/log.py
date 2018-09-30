"""
Utility functions for logging messages. Exports:

get_logger: Return a logger with a specified severity threshold.
"""

import logging
import sys
from . import constants




LOG_LEVEL = "log-level"


LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
}

