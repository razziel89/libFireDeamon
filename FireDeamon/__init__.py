import os

from .min import *
try:
    from .max import *
    FULL_SUPPORT = True
except ImportError as e:
    if os.environ.get("FD_REQ_FULL_SUPPORT", "0") == "1":
        raise e
    FULL_SUPPORT = False

