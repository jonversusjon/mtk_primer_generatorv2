from config.logging_config import logger  # Use the centralized logger
from typing import Union
from contextlib import contextmanager
import time
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
