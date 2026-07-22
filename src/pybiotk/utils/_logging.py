import logging
from typing import Union

from rich.console import Console
from rich.logging import RichHandler


def configure_logging(level: Union[int, str] = logging.INFO, rich: bool = False,
                      force: bool = False) -> None:
    """Configure standard logging, optionally using Rich for a CLI."""
    if rich:
        handlers = [
            RichHandler(
                console=Console(stderr=True, width=120),
                show_time=True,
                omit_repeated_times=False,
                show_level=True,
                markup=False,
                log_time_format="[%x %a %X]",
            )
        ]
        format_string = "%(message)s"
        date_format = None
    else:
        handlers = [logging.StreamHandler()]
        format_string = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        date_format = "%Y-%m-%d %H:%M:%S"

    logging.basicConfig(
        level=level,
        format=format_string,
        datefmt=date_format,
        handlers=handlers,
        force=force,
    )


def get_logger(name: str) -> logging.Logger:
    """Return a named logger after installing the default configuration."""
    configure_logging()
    return logging.getLogger(name)
