import logging
import time
import traceback
import json
import os
from contextlib import contextmanager
from config.logging_config import logger as base_logger  # Centralized logger
from pydantic import BaseModel
import numpy as np

class ModuleLoggerAdapter(logging.LoggerAdapter):
    def process(self, msg, kwargs):
        # Rely on Python's built-in module attribute rather than injecting our own.
        return msg, kwargs

class Logger:
    def __init__(self, name="MTKAdvanced", extra=None, enable_file_logging=True, log_dir="logs", custom_format=True):
        # Allow extra context info (e.g., {"module": __name__})
        self.extra = extra or {}
        # Create a child logger and wrap it with a LoggerAdapter for extra context.
        child_logger = base_logger.getChild(name)
        self.logger = ModuleLoggerAdapter(child_logger, self.extra)
        self.logger.logger.propagate = False  # Avoid duplicate logs
        self._setup_handlers(enable_file_logging, log_dir, custom_format)
        # Container for function timers.
        self.timers = {}

    def _setup_handlers(self, enable_file_logging, log_dir, custom_format):
        # If custom formatting is desired, remove any existing handlers and add our own.
        if custom_format:
            for handler in self.logger.logger.handlers[:]:
                self.logger.logger.removeHandler(handler)
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter(
                "%(asctime)s - %(levelname)s - %(message)s"
            )
            console_handler.setFormatter(console_formatter)
            self.logger.logger.addHandler(console_handler)
            self.logger.logger.setLevel(logging.DEBUG)
        # Optionally add a file handler.
        if enable_file_logging:
            os.makedirs(log_dir, exist_ok=True)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            file_handler = logging.FileHandler(
                f"{log_dir}/app_logger_{timestamp}.log", encoding="utf-8"
            )
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                "%(asctime)s | %(levelname)-8s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S"
            )
            file_handler.setFormatter(file_formatter)
            self.logger.logger.addHandler(file_handler)

    @contextmanager
    def debug_context(self, operation: str):
        """Context manager for debugging a block of code."""
        start_time = time.time()
        self.logger.debug(f"Starting: {operation}")
        try:
            yield
        except Exception as e:
            self.error(f"Error in {operation}: {str(e)}\n{traceback.format_exc()}", exc_info=True)
            raise
        finally:
            elapsed_time = time.time() - start_time
            self.logger.debug(f"Finished: {operation} (Time: {elapsed_time:.2f}s)")

    def log_function(self, func):
        """Decorator to log function entry, exit, parameters, and elapsed time."""
        def wrapper(*args, **kwargs):
            func_name = func.__name__
            self.logger.info(f"+{'-'*60}")
            self.logger.info(f"| STARTING: {func_name} with args={args} kwargs={kwargs}")
            self._start_timer(func_name)
            result = func(*args, **kwargs)
            elapsed = self._end_timer(func_name)
            self.logger.info(f"| COMPLETED: {func_name} in {elapsed:.4f} seconds with result={result}")
            self.logger.info(f"+{'-'*60}")
            return result
        return wrapper

    def _start_timer(self, name):
        self.timers[name] = time.time()

    def _end_timer(self, name):
        if name in self.timers:
            elapsed = time.time() - self.timers[name]
            del self.timers[name]
            return elapsed
        return None

    def log_step(self, step_name, message, data=None, level=logging.INFO):
        """Log a step or milestone in the code, ensuring JSON serialization of Pydantic models."""
        try:
            if isinstance(data, BaseModel):
                data_str = data.model_dump_json(indent=2)
            elif isinstance(data, list) and all(isinstance(i, BaseModel) for i in data):
                data_str = json.dumps([i.model_dump() for i in data], indent=2)
            elif isinstance(data, dict) and all(isinstance(v, BaseModel) for v in data.values()):
                data_str = json.dumps({k: v.model_dump() for k, v in data.items()}, indent=2)
            else:
                data_str = json.dumps(data, default=str, indent=2) if data else ""
        except (TypeError, ValueError) as e:
            data_str = f"[Failed to serialize data: {e}]"

        if data_str:
            log_message = f"{step_name} - {message}\nData: {data_str}"
        else:
            log_message = f"{step_name} - {message}"

        self.logger.log(level, log_message)

    @staticmethod
    def visualize_matrix(matrix, threshold=0):
        """
        Visualize a numpy matrix in ASCII format
        """
        if matrix.ndim <= 2:
            # 1D or 2D matrices
            if matrix.ndim == 1:
                matrix = matrix.reshape(1, -1)

            rows = []
            for row in matrix:
                row_str = " ".join(
                    ["#" if val > threshold else "." for val in row])
                rows.append(row_str)

            return "\n".join(rows)
        else:
            # For higher dimensions, show summary
            return f"Matrix shape: {matrix.shape}, non-zero: {np.count_nonzero(matrix)}"


    def validate(self, condition, message, data=None):
        """Log a validation result."""
        result = bool(condition)
        status = "PASS" if result else "FAIL"
        level = logging.INFO if result else logging.ERROR
        data_str = json.dumps(data, default=str, indent=2) if data else ""
        self.logger.log(level, f"VALIDATION {status}: {message} {data_str}")
        return result

    def debug(self, message, data=None):
        data_str = json.dumps(data, default=str, indent=2) if data else ""
        self.logger.debug(f"{message} {data_str}")

    # Added error method to allow logger.error calls
    def error(self, message, data=None, exc_info=False):
        data_str = json.dumps(data, default=str, indent=2) if data else ""
        self.logger.error(f"{message} {data_str}", exc_info=exc_info)

# Global instance for use throughout your app.
logger = Logger(log_dir="log_utils/logs")
