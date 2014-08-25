"""
Utility functions for logging messages. Exports:

get_logger: Return a logger with a specified severity threshold.
"""

import logging

LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
    }


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
    return logger
