import logging

from sargasso.utlis.log import LoggerManager


def _check_logger_level(logger, level):
    assert logger.getEffectiveLevel() == level


def test_get_logger_returns_debug_logger():
    logger = LoggerManager.get().init("debug").logger
    _check_logger_level(logger, logging.DEBUG)


def test_get_logger_returns_info_logger():
    logger = LoggerManager.get().init("info").logger
    _check_logger_level(logger, logging.INFO)


def test_get_logger_returns_warning_logger():
    logger = LoggerManager.get().init("warning").logger
    _check_logger_level(logger, logging.WARNING)


def test_get_logger_returns_error_logger():
    logger = LoggerManager.get().init("error").logger
    _check_logger_level(logger, logging.ERROR)


def test_get_logger_returns_critical_logger():
    logger = LoggerManager.get().init("critical").logger
    _check_logger_level(logger, logging.CRITICAL)
