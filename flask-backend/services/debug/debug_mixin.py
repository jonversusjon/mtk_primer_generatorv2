# debug_mixin.py
import logging
from functools import wraps

class DebugMixin:
    @staticmethod
    def debug_wrapper(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            if getattr(self, 'debugger', None):
                self.debugger.log_function_start(func.__name__, args, kwargs)
            result = func(self, *args, **kwargs)
            if getattr(self, 'debugger', None):
                self.debugger.log_function_end(func.__name__, result)
            return result
        return wrapper

    def log_step(self, step_name, message, data=None, level=logging.INFO):
        if getattr(self, 'debugger', None):
            self.debugger.log_step(step_name, message, data, level)

    def validate(self, condition, message, data=None):
        if getattr(self, 'debugger', None):
            return self.debugger.validate(condition, message, data)
        return condition
