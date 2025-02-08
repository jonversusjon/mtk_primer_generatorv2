# services/base.py
import logging
from typing import Optional, List
from contextlib import contextmanager
import time

class GoldenGateDesigner:
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
            formatter = logging.Formatter(
                '%(levelname)s - %(name)s - %(message)s'
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
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
        finally:
            duration = time.time() - start_time
            self.logger.debug(f"Completed {step_name} in {duration:.2f}s")

    def log_state(self, message: str):
        """Log current state with message."""
        self.logger.debug(f"{message} - Current state: {self.state}")