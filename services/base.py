# services/base.py
import logging
from typing import Optional, List
from contextlib import contextmanager
import time
from colorlog import ColoredFormatter
import traceback

class CustomFormatter(ColoredFormatter):
    def format(self, record):
        if record.levelname in ["DEBUG", "ERROR"]:
            self._fmt = "%(log_color)s%(levelname)s - %(name)s - %(filename)s:%(lineno)d - %(funcName)s - %(message)s"
        elif record.levelname == "INFO":
            self._fmt = "%(log_color)s%(message)s"  # No log level shown for INFO
        else:  # WARNING and others
            self._fmt = "%(log_color)s%(levelname)s - %(message)s"

        return super().format(record)
    
class PrimerDesignLogger:
    def __init__(self, verbose: bool = False, log_level: str = "DEBUG"):
        # Set up logging and state management
        self.logger = self._setup_logger(log_level)
        self.verbose = verbose
        self.state = {}  # Track execution state

    def _setup_logger(self, log_level: str) -> logging.Logger:
        """Creates a configured logger for the class."""
        logger = logging.getLogger(self.__class__.__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = CustomFormatter(
                "%(log_color)s%(levelname)s - %(message)s",
                log_colors={
                    'DEBUG': 'purple',
                    'INFO': 'white',
                    'WARNING': 'blue',
                    'ERROR': 'yellow',
                    'CRITICAL': 'red,bg_white'
                }
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(getattr(logging, log_level))
        return logger

    @contextmanager
    def debug_context(self, step_name: str):
        """Context manager for timing and logging operations."""
        self.logger.debug(f"Starting {step_name}")
        start_time = time.time()
        try:
            yield
        except Exception as e:
            tb = traceback.extract_tb(e.__traceback__)[-1]  # Get last frame of the traceback
            filename, line, func, text = tb
            self.logger.error(
                f"Error in {step_name} | Function: {func} | Line {line} | {e.__class__.__name__}: {str(e)}"
            )
            self.logger.debug("Traceback details:\n" + "".join(traceback.format_tb(e.__traceback__)))
            raise
        finally:
            duration = time.time() - start_time
            self.logger.debug(f"Completed {step_name} in {duration:.2f}s")


    def log_state(self, message: str):
        """Log current state with message."""
        self.logger.debug(f"{message} - Current state: {self.state}")