import logging
from piquant.log import get_logger

import sys


def _check_logger_level(logger, level):
    assert logger.getEffectiveLevel() == level


def test_get_logger_returns_debug_logger():
    logger = get_logger(sys.stdout, "debug")
    _check_logger_level(logger, logging.DEBUG)


def test_get_logger_returns_info_logger():
    logger = get_logger(sys.stdout, "info")
    _check_logger_level(logger, logging.INFO)


def test_get_logger_returns_warning_logger():
    logger = get_logger(sys.stdout, "warning")
    _check_logger_level(logger, logging.WARNING)


def test_get_logger_returns_error_logger():
    logger = get_logger(sys.stdout, "error")
    _check_logger_level(logger, logging.ERROR)


def test_get_logger_returns_critical_logger():
    logger = get_logger(sys.stdout, "critical")
    _check_logger_level(logger, logging.CRITICAL)
