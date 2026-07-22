import logging

from rich.logging import RichHandler

from pybiotk.utils import configure_logging, get_logger
from pybiotk.utils import logging as compatibility_logging


def test_logging_compatibility_export():
    assert compatibility_logging is logging


def test_get_logger_uses_standard_default_configuration():
    root_logger = logging.getLogger()
    original_handlers = root_logger.handlers[:]
    original_level = root_logger.level
    try:
        root_logger.handlers = []
        logger = get_logger("test.module")

        handler = root_logger.handlers[0]
        assert type(handler) is logging.StreamHandler
        assert logger.name == "test.module"
        assert handler.formatter._fmt == "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        assert handler.formatter.datefmt == "%Y-%m-%d %H:%M:%S"
    finally:
        root_logger.handlers = original_handlers
        root_logger.setLevel(original_level)


def test_configure_logging_uses_rich_stderr_console():
    root_logger = logging.getLogger()
    original_handlers = root_logger.handlers[:]
    original_level = root_logger.level
    try:
        configure_logging(rich=True, force=True)

        handler = root_logger.handlers[0]
        assert isinstance(handler, RichHandler)
        assert handler.console.stderr is True
        assert handler.console.width == 120
    finally:
        root_logger.handlers = original_handlers
        root_logger.setLevel(original_level)
