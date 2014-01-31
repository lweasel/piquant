import logging

LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
    }


def getLogger(stream, level):
    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s: %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler(stream)
    handler.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    logger.setLevel(LEVELS[level])
    logger.addHandler(handler)
    return logger
