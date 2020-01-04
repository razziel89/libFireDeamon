from .cpp_min import *
try:
    from .cpp_max import *
    MIN_SUPPORT = False
except ImportError:
    MIN_SUPPORT = True
