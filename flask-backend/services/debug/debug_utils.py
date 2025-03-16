import logging
import json
import time
import os
import numpy as np

from config.logging_config import logger  # Use the centralized logger
from contextlib import contextmanager
import traceback

@contextmanager
def debug_context(operation: str):
    """Context manager for debugging operations."""
    start_time = time.time()
    logger.debug(f"Starting: {operation}")
    try:
        yield
    except Exception as e:
        logger.error(f"Error in {operation}: {str(e)}\n{traceback.format_exc()}")
        raise
    finally:
        elapsed_time = time.time() - start_time
        logger.debug(f"Finished: {operation} (Time: {elapsed_time:.2f}s)")


class MutationDebugger:
    """
    A comprehensive debugging utility for the MutationOptimizer class that
    integrates with existing logging configuration.
    """

    def __init__(self,
                 parent_logger=None,
                 log_dir="logs",
                 enable_file_logging=True,
                 use_custom_format=True):
        """
        Initialize the debugger with customizable logging options.

        Args:
            parent_logger: The parent logger to use (if None, creates a new one)
            log_dir (str): Directory to store log files
            enable_file_logging (bool): Whether to log to files in addition to console
            use_custom_format (bool): Whether to override application's log format for cleaner output
        """
        if parent_logger:
            # Use the provided parent logger
            self.logger = parent_logger.getChild("Debugger")
        else:
            # Create new logger
            self.logger = logging.getLogger("MutationDebugger")

        # Prevent propagation to avoid duplicate logs
        self.logger.propagate = False

        # Store all existing handlers to restore later if needed
        self._existing_handlers = list(self.logger.handlers)

        # Remove existing handlers if using custom format
        if use_custom_format:
            for handler in self.logger.handlers[:]:
                self.logger.removeHandler(handler)

            # Add console handler with clean format
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter("%(message)s")
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)

            # Ensure the logger level is set to DEBUG
            self.logger.setLevel(logging.DEBUG)

        # Create file handler if enabled
        if enable_file_logging:
            os.makedirs(log_dir, exist_ok=True)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            file_handler = logging.FileHandler(
                f"{log_dir}/mutation_debug_{timestamp}.log",
                encoding='utf-8'
            )
            file_handler.setLevel(logging.DEBUG)

            # Use detailed format for file logs
            file_formatter = logging.Formatter(
                "%(asctime)s | %(levelname)-8s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S"
            )
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)

        self.step_counts = {}
        self.timers = {}
        self.validation_results = {}
        self.current_function = None

    def __del__(self):
        """Clean up by restoring original handlers if needed"""
        try:
            # Only attempt cleanup if the logger still exists
            if hasattr(self, 'logger') and hasattr(self, '_existing_handlers'):
                # Remove any handlers we added
                for handler in self.logger.handlers[:]:
                    self.logger.removeHandler(handler)

                # Restore original handlers
                for handler in self._existing_handlers:
                    self.logger.addHandler(handler)
        except Exception as e:  # Changed 'catch' to 'except' and added 'Exception as'
            print(f"Error during cleanup: {e}")
            
    def start_timer(self, name):
        """Start a timer for performance tracking"""
        self.timers[name] = time.time()

    def end_timer(self, name):
        """End a timer and return the elapsed time"""
        if name in self.timers:
            elapsed = time.time() - self.timers[name]
            del self.timers[name]
            return elapsed
        return None

    def log_function_start(self, func_name, args=None, kwargs=None):
        """Log the start of a function with its arguments"""
        self.current_function = func_name
        self.start_timer(func_name)

        # Format args and kwargs for pretty printing
        args_str = ""
        if args:
            args_str = ", ".join([self._format_arg(arg) for arg in args])

        kwargs_str = ""
        if kwargs:
            kwargs_str = ", ".join(
                [f"{k}={self._format_arg(v)}" for k, v in kwargs.items()])

        params = ", ".join(filter(None, [args_str, kwargs_str]))

        self.logger.info(f"+{'-' * 60}")
        self.logger.info(f"| STARTING: {func_name}({params})")
        self.logger.info(f"+{'-' * 60}")

    def log_function_end(self, func_name, result=None):
        """Log the end of a function with its result summary"""
        elapsed = self.end_timer(func_name)

        self.logger.info(f"+{'-' * 60}")
        if result is not None:
            result_summary = self._summarize_result(result)
            self.logger.info(f"| RESULT: {result_summary}")

        self.logger.info(f"| COMPLETED: {func_name} in {elapsed:.4f} seconds")
        self.logger.info(f"+{'-' * 60}")
        self.current_function = None

    def log_step(self, step_name, message, data=None, level=logging.INFO):
        """Log a step within a function with optional data"""
        if self.current_function not in self.step_counts:
            self.step_counts[self.current_function] = 0

        self.step_counts[self.current_function] += 1
        step_num = self.step_counts[self.current_function]

        prefix = f"| STEP {step_num}: {step_name}"
        self.logger.log(level, f"{prefix} - {message}")

        if data is not None:
            data_str = self._format_data(data)
            for line in data_str.split('\n'):
                self.logger.log(level, f"|   {line}")

    def validate(self, condition, message, data=None):
        """Validate a condition and log the result"""
        result = bool(condition)
        status = "PASS" if result else "FAIL"

        if self.current_function not in self.validation_results:
            self.validation_results[self.current_function] = []

        self.validation_results[self.current_function].append(result)

        log_level = logging.INFO if result else logging.ERROR
        self.logger.log(log_level, f"| VALIDATION {status}: {message}")

        if not result and data is not None:
            data_str = self._format_data(data)
            for line in data_str.split('\n'):
                self.logger.log(log_level, f"|   {line}")

        return result

    def log_warning(self, message, data=None):
        """Log a warning with optional data"""
        self.logger.warning(f"| WARNING: {message}")
        if data is not None:
            data_str = self._format_data(data)
            for line in data_str.split('\n'):
                self.logger.warning(f"|   {line}")

    def log_error(self, message, data=None):
        """Log an error with optional data"""
        self.logger.error(f"| ERROR: {message}")
        if data is not None:
            data_str = self._format_data(data)
            for line in data_str.split('\n'):
                self.logger.error(f"|   {line}")

    def _format_arg(self, arg):
        """Format an argument for logging"""
        if isinstance(arg, str):
            if len(arg) > 50:
                return f"'{arg[:47]}...'"
            return f"'{arg}'"
        elif isinstance(arg, (list, tuple)):
            if len(arg) > 3:
                return f"{type(arg).__name__} of {len(arg)} items"
            return str(arg)
        elif isinstance(arg, dict):
            if len(arg) > 3:
                return f"dict with {len(arg)} keys"
            return str(arg)
        elif isinstance(arg, np.ndarray):
            return f"ndarray shape={arg.shape}"
        else:
            return str(arg)

    def _summarize_result(self, result):
        """Summarize a result for logging"""
        if isinstance(result, tuple) and len(result) == 2:
            # Handle case where there are two return values
            result1, result2 = result
            summary1 = self._format_arg(result1)
            summary2 = self._format_arg(result2)
            return f"({summary1}, {summary2})"
        else:
            return self._format_arg(result)

    def _format_data(self, data):
        """Format data for logging"""
        if isinstance(data, dict):
            return json.dumps(data, indent=2, default=str)
        elif isinstance(data, (list, tuple)):
            if len(data) > 5:
                items = [str(item) for item in data[:5]]
                return f"[{', '.join(items)}, ... ({len(data)-5} more)]"
            return str(data)
        elif isinstance(data, np.ndarray):
            if data.size > 25:
                return f"ndarray shape={data.shape}, dtype={data.dtype}, non-zero={np.count_nonzero(data)}"
            return str(data)
        else:
            return str(data)

    def summarize_validations(self):
        """Summarize all validation results"""
        self.logger.info("-" * 70)
        self.logger.info("VALIDATION SUMMARY")
        self.logger.info("-" * 70)

        all_pass = True

        for func_name, results in self.validation_results.items():
            pass_count = sum(results)
            total_count = len(results)

            if pass_count == total_count:
                status = "ALL PASS"
            else:
                status = f"{total_count - pass_count} FAIL"
                all_pass = False

            self.logger.info(
                f"{func_name}: {status} ({pass_count}/{total_count})")

        if all_pass:
            self.logger.info("ALL VALIDATIONS PASSED")
        else:
            self.logger.error("SOME VALIDATIONS FAILED")

        self.logger.info("-" * 70)

        return all_pass


def debug_function(func):
    """
    Decorator for automatically logging function entry/exit
    This is a marker decorator that will be replaced at runtime
    """
    return func


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


def visualize_overhang_compatibility(overhangs, compatibility_table, utils):
    """
    Visualize the compatibility of a set of overhangs
    """
    n = len(overhangs)
    result = []

    # Header
    header = "     | " + " ".join([f"{i+1:4}" for i in range(n)])
    result.append(header)
    result.append("-" * len(header))

    # Compatibility matrix
    for i in range(n):
        row = [f"{i+1:4} | "]
        for j in range(n):
            if i == j:
                row.append("  .  ")
                continue

            idx1 = utils.seq_to_index(overhangs[i])
            idx2 = utils.seq_to_index(overhangs[j])

            if compatibility_table[idx1][idx2] == 1:
                row.append("  +  ")
            else:
                row.append("  -  ")

        result.append("".join(row))

    # Add overhang sequences for reference
    result.append("-" * len(header))
    for i, overhang in enumerate(overhangs):
        result.append(f"{i+1:4} : {overhang}")

    return "\n".join(result)
